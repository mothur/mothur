/*
 *  aligncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 *	This version of nast does everything I think that the greengenes nast server does and then some.  I have added the 
 *	feature of allowing users to define their database, kmer size for searching, alignment penalty values and alignment 
 *	method.  This latter feature is perhaps most significant.  nastPlus enables a user to use either a Needleman-Wunsch 
 *	(non-affine gap penalty) or Gotoh (affine gap penalty) pairwise alignment algorithm.  This is significant because it
 *	allows for a global alignment and not the local alignment provided by bLAst.  Furthermore, it has the potential to
 *	provide a better alignment because of the banding method employed by blast (I'm not sure about this).
 *
 */

#include "aligncommand.h"
#include "sequence.hpp"

#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "kmerdb.hpp"
#include "suffixdb.hpp"
#include "blastdb.hpp"

#include "nast.hpp"
#include "nastreport.hpp"


//**********************************************************************************************************************

AlignCommand::AlignCommand(string option){
	try {
		//		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"template","candidate","search","ksize","align","match","mismatch","gapopen","gapextend", "processors"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			templateFileName = validParameter.validFile(parameters, "template", true);
			if (templateFileName == "not found") { 
				mothurOut("template is a required parameter for the align.seqs command."); 
				mothurOutEndLine();
				abort = true; 
			}
			else if (templateFileName == "not open") { abort = true; }	
			
			candidateFileName = validParameter.validFile(parameters, "candidate", true);
			if (candidateFileName == "not found") { 
				mothurOut("candidate is a required parameter for the align.seqs command."); 
				mothurOutEndLine();
				abort = true; 
			}
			else if (candidateFileName == "not open") { abort = true; }	
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "ksize", false);		if (temp == "not found"){	temp = "8";				}
			convert(temp, kmerSize); 
			
			temp = validParameter.validFile(parameters, "match", false);		if (temp == "not found"){	temp = "1.0";			}
			convert(temp, match);  
			
			temp = validParameter.validFile(parameters, "mismatch", false);		if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, misMatch);  
			
			temp = validParameter.validFile(parameters, "gapopen", false);		if (temp == "not found"){	temp = "-2.0";			}
			convert(temp, gapOpen);  
			
			temp = validParameter.validFile(parameters, "gapextend", false);	if (temp == "not found"){	temp = "-1.0";			}
			convert(temp, gapExtend); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			search = validParameter.validFile(parameters, "search", false);		if (search == "not found"){	search = "kmer";		}
			
			align = validParameter.validFile(parameters, "align", false);		if (align == "not found"){	align = "needleman";	}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "AlignCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

AlignCommand::~AlignCommand(){	

	if (abort == false) {
		for (int i = 0; i < lines.size(); i++) {  delete lines[i];  }  lines.clear();
		delete templateDB;
		delete alignment;
	}
}

//**********************************************************************************************************************

void AlignCommand::help(){
	try {
		mothurOut("The align.seqs command reads a file containing sequences and creates an alignment file and a report file.\n");
		mothurOut("The align.seqs command parameters are template, candidate, search, ksize, align, match, mismatch, gapopen and gapextend.\n");
		mothurOut("The template and candidate parameters are required.\n");
		mothurOut("The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer.\n");
		mothurOut("The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman.\n");
		mothurOut("The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 7.\n");
		mothurOut("The match parameter allows you to specify the bonus for having the same base. The default is 1.0.\n");
		mothurOut("The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0.\n");
		mothurOut("The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -1.0.\n");
		mothurOut("The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -2.0.\n");
		mothurOut("The align.seqs command should be in the following format: \n");
		mothurOut("align.seqs(template=yourTemplateFile, candidate=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty) \n");
		mothurOut("Example align.seqs(candidate=candidate.fasta, template=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)\n");
		mothurOut("Note: No spaces between parameter labels (i.e. candidate), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int AlignCommand::execute(){
	try {
		if (abort == true) {	return 0;	}
		
		if(search == "kmer")			{	templateDB = new KmerDB(templateFileName, kmerSize);	}
		else if(search == "suffix")		{	templateDB = new SuffixDB(templateFileName);			}
		else if(search == "blast")		{	templateDB = new BlastDB(templateFileName, gapOpen, gapExtend, match, misMatch);	}
		else {
			mothurOut(search + " is not a valid search option. I will run the command using kmer, ksize=8.");
			mothurOutEndLine();
			kmerSize = 8;
			
			templateDB = new KmerDB(templateFileName, kmerSize);
		}
		
		int longestBase = templateDB->getLongestBase();
		
		if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, longestBase);			}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);				}
		else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
		else if(align == "noalign")		{	alignment = new NoAlign();													}
		else {
			mothurOut(align + " is not a valid alignment option. I will run the command using needleman.");
			mothurOutEndLine();
			alignment = new NeedlemanOverlap(gapOpen, match, misMatch, longestBase);
		}
		
		string alignFileName = candidateFileName.substr(0,candidateFileName.find_last_of(".")+1) + "align";
		string reportFileName = candidateFileName.substr(0,candidateFileName.find_last_of(".")+1) + "align.report";
		
		int numFastaSeqs = 0;
		int start = time(NULL);
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1){
			ifstream inFASTA;
			openInputFile(candidateFileName, inFASTA);
			numFastaSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
			inFASTA.close();
			
			lines.push_back(new linePair(0, numFastaSeqs));
		
			driver(lines[0], alignFileName, reportFileName);
			
		}
		else{
			vector<int> positions;
			processIDS.resize(0);
			
			ifstream inFASTA;
			openInputFile(candidateFileName, inFASTA);
			
			string input;
			while(!inFASTA.eof()){
				input = getline(inFASTA);
				if (input.length() != 0) {
					if(input[0] == '>'){	int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
				}
			}
			inFASTA.close();
			
			numFastaSeqs = positions.size();
			
			int numSeqsPerProcessor = numFastaSeqs / processors;
			
			for (int i = 0; i < processors; i++) {
				int startPos = positions[ i * numSeqsPerProcessor ];
				if(i == processors - 1){
					numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;
				}
				lines.push_back(new linePair(startPos, numSeqsPerProcessor));
			}
			createProcesses(alignFileName, reportFileName); 
			
			rename((alignFileName + toString(processIDS[0]) + ".temp").c_str(), alignFileName.c_str());
			rename((reportFileName + toString(processIDS[0]) + ".temp").c_str(), reportFileName.c_str());
			
			for(int i=1;i<processors;i++){
				appendAlignFiles((alignFileName + toString(processIDS[i]) + ".temp"), alignFileName);
				remove((alignFileName + toString(processIDS[i]) + ".temp").c_str());
				
				appendReportFiles((reportFileName + toString(processIDS[i]) + ".temp"), reportFileName);
				remove((reportFileName + toString(processIDS[i]) + ".temp").c_str());
			}
			
		}
#else
		ifstream inFASTA;
		openInputFile(candidateFileName, inFASTA);
		numFastaSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
		inFASTA.close();
		
		lines.push_back(new linePair(0, numFastaSeqs));
		
		driver(lines[0], alignFileName, reportFileName);
#endif
		
		mothurOut("It took " + toString(time(NULL) - start) + " secs to align " + toString(numFastaSeqs) + " sequences.");
		mothurOutEndLine();
		mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************

int AlignCommand::driver(linePair* line, string alignFName, string reportFName){
	try {
		ofstream alignmentFile;
		openOutputFile(alignFName, alignmentFile);
		NastReport report(reportFName);
		
		ifstream inFASTA;
		openInputFile(candidateFileName, inFASTA);

		inFASTA.seekg(line->start);

		for(int i=0;i<line->numSeqs;i++){
			
			Sequence* candidateSeq = new Sequence(inFASTA);
			report.setCandidate(candidateSeq);
	
			Sequence temp = templateDB->findClosestSequence(candidateSeq);
			Sequence* templateSeq = &temp;
			
			report.setTemplate(templateSeq);
			report.setSearchParameters(search, templateDB->getSearchScore());
			
			Nast nast(alignment, candidateSeq, templateSeq);

			report.setAlignmentParameters(align, alignment);

			report.setNastParameters(nast);
			
			alignmentFile << '>' << candidateSeq->getName() << '\n' << candidateSeq->getAligned() << endl;
			
			report.print();
			
			delete candidateSeq;		
		}
		
		alignmentFile.close();
		inFASTA.close();
		
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "driver");
		exit(1);
	}
}

/**************************************************************************************************/

void AlignCommand::createProcesses(string alignFileName, string reportFileName) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		//		processIDS.resize(0);
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driver(lines[process], alignFileName + toString(getpid()) + ".temp", reportFileName + toString(getpid()) + ".temp");
				exit(0);
			}else { mothurOut("unable to spawn the necessary processes."); mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
#endif		
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "createProcesses");
		exit(1);
	}
}

/**************************************************************************************************/

void AlignCommand::appendAlignFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		openOutputFileAppend(filename, output);
		openInputFile(temp, input);
		
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
		}
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "appendAlignFiles");
		exit(1);
	}
}

/**************************************************************************************************/

void AlignCommand::appendReportFiles(string temp, string filename) {
	try{
		
		ofstream output;
		ifstream input;
		openOutputFileAppend(filename, output);
		openInputFile(temp, input);

		while (!input.eof())	{	char c = input.get(); if (c == 10 || c == 13){	break;	}	} // get header line
				
		while(char c = input.get()){
			if(input.eof())		{	break;			}
			else				{	output << c;	}
		}
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		errorOut(e, "AlignCommand", "appendReportFiles");
		exit(1);
	}
}

//**********************************************************************************************************************
