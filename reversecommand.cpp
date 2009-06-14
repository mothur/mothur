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

ReverseSeqsCommand::ReverseSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fasta = validParameter.validFile(parameters, "fasta", true);
			if (fasta == "not open") { abort = true; }
			else if (fasta == "not found") { fasta = ""; cout << "fasta is a required parameter for the reverse.seqs command." << endl; abort = true;  }	
			
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReverseSeqsCommand class Function ReverseSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReverseSeqsCommand class function ReverseSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

void ReverseSeqsCommand::help(){
	try {
		cout << "The reverse.seqs command reads a fastafile and ...." << "\n";
		cout << "The reverse.seqs command parameter is fasta and it is required." << "\n";
		cout << "The reverse.seqs command should be in the following format: " << "\n";
		cout << "reverse.seqs(fasta=yourFastaFile) " << "\n";	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReverseSeqsCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReverseSeqsCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//***************************************************************************************************************

ReverseSeqsCommand::~ReverseSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************


int ReverseSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		ifstream inFASTA;
		openInputFile(fasta, inFASTA);
		
		ofstream outFASTA;
		string reverseFile = getRootName(fasta) + "rc" + getExtension(fasta);
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
		cout << "Standard Error: " << e.what() << " has occurred in the ReverseSeqsCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReverseSeqsCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//***************************************************************************************************************
