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
			if (templateFileName == "not found") { cout << "template is a required parameter for the align.seqs command." << endl; abort = true; }
			else if (templateFileName == "not open") { abort = true; }	
			
			candidateFileName = validParameter.validFile(parameters, "candidate", true);
			if (candidateFileName == "not found") { cout << "candidate is a required parameter for the align.seqs command." << endl; abort = true; }
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
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function AlignCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function AlignCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

AlignCommand::~AlignCommand(){			
	delete templateDB;
	delete alignment;
}

//**********************************************************************************************************************

void AlignCommand::help(){
	try {
		cout << "The align.seqs command reads a file containing sequences and creates an alignment file and a report file." << "\n";
		cout << "The align.seqs command parameters are template, candidate, search, ksize, align, match, mismatch, gapopen and gapextend.  " << "\n";
		cout << "The template and candidate parameters are required." << "\n";
		cout << "The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer." << "\n";
		cout << "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman." << "\n";
		cout << "The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 7." << "\n";
		cout << "The match parameter allows you to specify the bonus for having the same base. The default is 1.0." << "\n";
		cout << "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0." << "\n";
		cout << "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -1.0." << "\n";
		cout << "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -2.0." << "\n";
		cout << "The align.seqs command should be in the following format: " << "\n";
		cout << "align.seqs(template=yourTemplateFile, candidate=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty) " << "\n";
		cout << "Example align.seqs(candidate=candidate.fasta, template=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)" << "\n";
		cout << "Note: No spaces between parameter labels (i.e. candidate), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
			cout << search << " is not a valid search option. I will run the command using kmer, ksize=8." << endl; kmerSize = 8;
			templateDB = new KmerDB(templateFileName, kmerSize);
		}
		
		if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, 3000);	}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 3000);			}
		else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
		else if(align == "noalign")		{	alignment = new NoAlign();													}
		else {
			cout << align << " is not a valid alignment option. I will run the command using needleman." << endl;
			alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 3000);
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
			
			while(!inFASTA.eof()){
				char c = inFASTA.get();
				if(c == '>'){	positions.push_back(inFASTA.tellg());	}
				while (!inFASTA.eof())	{	c = inFASTA.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
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
		
		cout << "It took " << time(NULL) - start << " secs to align " << numFastaSeqs << " sequences" << endl;
		cout << endl;
		
		
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
			
			Sequence* templateSeq = templateDB->findClosestSequence(candidateSeq);
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
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function driver. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function driver. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
			}else { cout << "unable to spawn the necessary processes." << endl; exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
#endif		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the AlignCommand class Function createProcesses. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the AlignCommand class function createProcesses. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/**************************************************************************************************/

void AlignCommand::appendAlignFiles(string temp, string filename) {
	try{
		
		//open output file in append mode
		ofstream output;
		openOutputFileAppend(filename, output);
		
		//open temp file for reading
		ifstream input;
		openInputFile(temp, input);
		
		string line;
		//read input file and write to output file
		while(input.eof() != true) {
			getline(input, line); //getline removes the newline char
			if (line != "") {
				output << line << endl;   // Appending back newline char 
			}
		}	
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function appendFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function appendFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/**************************************************************************************************/

void AlignCommand::appendReportFiles(string temp, string filename) {
	try{
		
		//open output file in append mode
		ofstream output;
		openOutputFileAppend(filename, output);
		
		//open temp file for reading
		ifstream input;
		openInputFile(temp, input);
		while (!input.eof())	{	char c = input.get(); if (c == 10 || c == 13){	break;	}	} // get header line
		
		string line;
		//read input file and write to output file
		while(input.eof() != true) {
			getline(input, line); //getline removes the newline char
			if (line != "") {
				output << line << endl;   // Appending back newline char 
			}
		}	
		
		input.close();
		output.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the DistanceCommand class Function appendFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the DistanceCommand class function appendFiles. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
