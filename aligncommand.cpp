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
 *	to compile type:
 *		make
 *
 *	for basic instructions on how to run nastPlus type:
 *		./nastPlus
 */

#include "aligncommand.h"
#include "sequence.hpp"

#include "alignment.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "database.hpp"
#include "kmerdb.hpp"
#include "suffixdb.hpp"
#include "blastdb.hpp"

#include "nast.hpp"
#include "nastreport.hpp"


//**********************************************************************************************************************
AlignCommand::AlignCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta","candidate","search","ksize","align","match","mismatch","gapopen","gapextend"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			parser = new OptionParser();

			parser->parse(option, parameters);   	delete parser; 
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter->isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			templateFileName = validParameter->validFile(parameters, "fasta", true);
			if (templateFileName == "not found") { cout << "fasta is a required parameter for the align.seqs command." << endl; abort = true; }
			else if (templateFileName == "not open") { abort = true; }	
			else { globaldata->setFastaFile(templateFileName); }
		
			candidateFileName = validParameter->validFile(parameters, "candidate", true);
			if (candidateFileName == "not found") { cout << "candidate is a required parameter for the align.seqs command." << endl; abort = true; }
			else if (candidateFileName == "not open") { abort = true; }	
			else { 
				globaldata->setCandidateFile(candidateFileName);
				openInputFile(candidateFileName, in);		
			}
			
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter->validFile(parameters, "ksize", false);			if (temp == "not found") { temp = "8"; }
			convert(temp, kmerSize); 
		
			temp = validParameter->validFile(parameters, "match", false);			if (temp == "not found") { temp = "1.0"; }
			convert(temp, match);  

			temp = validParameter->validFile(parameters, "mismatch", false);			if (temp == "not found") { temp = "-1.0"; }
			convert(temp, misMatch);  

			temp = validParameter->validFile(parameters, "gapopen", false);			if (temp == "not found") { temp = "-1.0"; }
			convert(temp, gapOpen);  

			temp = validParameter->validFile(parameters, "gapextend", false);		if (temp == "not found") { temp = "-2.0"; }
			convert(temp, gapExtend); 
		
			search = validParameter->validFile(parameters, "search", false);			if (search == "not found")	{ search = "kmer";		}
			align = validParameter->validFile(parameters, "align", false);			if (align == "not found")	{ align = "needleman";	}
			
			delete validParameter;
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
	
}
//**********************************************************************************************************************

void AlignCommand::help(){
	try {
		cout << "The align.seqs command reads a file containing sequences and creates an alignment file and a report file." << "\n";
		cout << "The align.seqs command parameters are fasta, candidate, search, ksize, align, match, mismatch, gapopen and gapextend.  " << "\n";
		cout << "The fasta and candidate parameters are required." << "\n";
		cout << "The search parameter allows you to specify the method to find most similar template.  Your options are: suffix, kmer and blast. The default is kmer." << "\n";
		cout << "The align parameter allows you to specify the alignment method to use.  Your options are: gotoh, needleman, blast and noalign. The default is needleman." << "\n";
		cout << "The ksize parameter allows you to specify the kmer size for finding most similar template to candidate.  The default is 7." << "\n";
		cout << "The match parameter allows you to specify the bonus for having the same base. The default is 1.0." << "\n";
		cout << "The mistmatch parameter allows you to specify the penalty for having different bases.  The default is -1.0." << "\n";
		cout << "The gapopen parameter allows you to specify the penalty for opening a gap in an alignment. The default is -1.0." << "\n";
		cout << "The gapextend parameter allows you to specify the penalty for extending a gap in an alignment.  The default is -2.0." << "\n";
		cout << "The align.seqs command should be in the following format: " << "\n";
		cout << "align.seqs(fasta=yourTemplateFile, candidate=yourCandidateFile, align=yourAlignmentMethod, search=yourSearchmethod, ksize=yourKmerSize, match=yourMatchBonus, mismatch=yourMismatchpenalty, gapopen=yourGapopenPenalty, gapextend=yourGapExtendPenalty) " << "\n";
		cout << "Example align.seqs(candidate=candidate.fasta, fasta=core.filtered, align=kmer, search=gotoh, ksize=8, match=2.0, mismatch=3.0, gapopen=-2.0, gapextend=-1.0)" << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
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
	
		srand( (unsigned)time( NULL ) );  //needed to assign names to temporary files
		
		Database* templateDB;
		if(search == "kmer")			{	templateDB = new KmerDB(templateFileName, kmerSize);	}
		else if(search == "suffix")		{	templateDB = new SuffixDB(templateFileName);			}
		else if(search == "blast")		{	templateDB = new BlastDB(templateFileName, gapOpen, gapExtend, match, misMatch);	}
		else {
			cout << search << " is not a valid search option. I will run the command using kmer, ksize=8." << endl; kmerSize = 8;
			templateDB = new KmerDB(templateFileName, kmerSize);
		}
	
		Alignment* alignment;
		if(align == "gotoh")			{	alignment = new GotohOverlap(gapOpen, gapExtend, match, misMatch, 3000);	}
		else if(align == "needleman")	{	alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 3000);			}
		else if(align == "blast")		{	alignment = new BlastAlignment(gapOpen, gapExtend, match, misMatch);		}
		else if(align == "noalign")		{	alignment = new NoAlign();													}
		else {
			cout << align << " is not a valid alignment option. I will run the command using needleman." << endl;
			alignment = new NeedlemanOverlap(gapOpen, match, misMatch, 3000);
		}
				
		int numFastaSeqs=count(istreambuf_iterator<char>(in),istreambuf_iterator<char>(), '>');
		in.seekg(0);
	
		string candidateAligngmentFName = candidateFileName.substr(0,candidateFileName.find_last_of(".")+1) + "align";
		ofstream candidateAlignmentFile;
		openOutputFile(candidateAligngmentFName, candidateAlignmentFile);

		string candidateReportFName = candidateFileName.substr(0,candidateFileName.find_last_of(".")+1) + "align.report";
		NastReport report(candidateReportFName);

		cout << "We are going to align the " << numFastaSeqs << " sequences in " << candidateFileName << "..." << endl;
		cout.flush();
	
		int start = time(NULL);

		for(int i=0;i<numFastaSeqs;i++){

			Sequence* candidateSeq = new Sequence(in);
			report.setCandidate(candidateSeq);

			Sequence* templateSeq = templateDB->findClosestSequence(candidateSeq);
			report.setTemplate(templateSeq);
			report.setSearchParameters(search, templateDB->getSearchScore());
			
				
			Nast nast(alignment, candidateSeq, templateSeq);
			report.setAlignmentParameters(align, alignment);
			report.setNastParameters(nast);

			candidateAlignmentFile << '>' << candidateSeq->getName() << '\n' << candidateSeq->getAligned() << endl;
			candidateAlignmentFile.flush();

			if((i+1) % 100 == 0){
				cout << "It has taken " << time(NULL) - start << " secs to align " << i+1 << " sequences" << endl;
			}
			report.print();
		
			delete candidateSeq;		
		}
		
		cout << "It took " << time(NULL) - start << " secs to align " << numFastaSeqs << " sequences" << endl;
		cout << endl;

		delete templateDB;
		delete alignment;

				
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