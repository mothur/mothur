/*
 *  reversecommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "reversecommand.h"
#include "sequence.hpp"


//***************************************************************************************************************

ReverseSeqsCommand::ReverseSeqsCommand(){
	try {
		globaldata = GlobalData::getInstance();
		if(globaldata->getFastaFile() == "")		{	cout << "you need to at least enter a fasta file name" << endl;	}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SeqCoordCommand class Function SeqCoordCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SeqCoordCommand class function SeqCoordCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//***************************************************************************************************************

ReverseSeqsCommand::~ReverseSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************


int ReverseSeqsCommand::execute(){
	try{
		
		ifstream inFASTA;
		openInputFile(globaldata->getFastaFile(), inFASTA);
		
		ofstream outFASTA;
		string reverseFile = getRootName(globaldata->getFastaFile()) + "rc" + getExtension(globaldata->getFastaFile());
		openOutputFile(reverseFile, outFASTA);
		
		while(!inFASTA.eof()){
			Sequence currSeq(inFASTA);
			currSeq.reverseComplement();
			currSeq.printSequence(outFASTA);
		}
		inFASTA.close();
		outFASTA.close();
		
		return 0;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the FilterSeqsCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the FilterSeqsCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//***************************************************************************************************************
