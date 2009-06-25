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
			else if (fasta == "not found") { fasta = ""; mothurOut("fasta is a required parameter for the reverse.seqs command."); mothurOutEndLine(); abort = true;  }	
			
		}
	}
	catch(exception& e) {
		errorOut(e, "ReverseSeqsCommand", "ReverseSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ReverseSeqsCommand::help(){
	try {
		mothurOut("The reverse.seqs command reads a fastafile and ....\n");
		mothurOut("The reverse.seqs command parameter is fasta and it is required.\n");
		mothurOut("The reverse.seqs command should be in the following format: \n");
		mothurOut("reverse.seqs(fasta=yourFastaFile) \n");	
	}
	catch(exception& e) {
		errorOut(e, "ReverseSeqsCommand", "help");
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
		errorOut(e, "ReverseSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************
