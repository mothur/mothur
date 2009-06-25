/*
 *  deconvolute.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "deconvolutecommand.h"

/**************************************************************************************/
DeconvoluteCommand::DeconvoluteCommand(string option) {	
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "name"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			inFastaName = validParameter.validFile(parameters, "fasta", true);
			if (inFastaName == "not open") { abort = true; }
			else if (inFastaName == "not found") { inFastaName = ""; mothurOut("fasta is a required parameter for the unique.seqs command."); mothurOutEndLine(); abort = true;  }	
			
			oldNameMapFName = validParameter.validFile(parameters, "name", true);
			if (oldNameMapFName == "not open") { abort = true; }
			else if (oldNameMapFName == "not found"){	oldNameMapFName = "";	}
		}

	}
	catch(exception& e) {
		errorOut(e, "DeconvoluteCommand", "DeconvoluteCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void DeconvoluteCommand::help(){
	try {
		mothurOut("The unique.seqs command reads a fastafile and creates a namesfile.\n");
		mothurOut("It creates a file where the first column is the groupname and the second column is a list of sequence names who have the same sequence. \n");
		mothurOut("If the sequence is unique the second column will just contain its name. \n");
		mothurOut("The unique.seqs command parameter is fasta and it is required.\n");
		mothurOut("The unique.seqs command should be in the following format: \n");
		mothurOut("unique.seqs(fasta=yourFastaFile) \n");	
	}
	catch(exception& e) {
		errorOut(e, "DeconvoluteCommand", "help");
		exit(1);
	}
}

/**************************************************************************************/
int DeconvoluteCommand::execute() {	
	try {
		
		if (abort == true) { return 0; }

		//prepare filenames and open files
		string outNameFile = (getRootName(inFastaName) + "names");
		string outFastaFile = (getRootName(inFastaName) + "unique" + getExtension(inFastaName));
		
		FastaMap fastamap;
	
		if(oldNameMapFName == "")	{	fastamap.readFastaFile(inFastaName);					}
		else						{	fastamap.readFastaFile(inFastaName, oldNameMapFName);	}
		
		fastamap.printCondensedFasta(outFastaFile);
		fastamap.printNamesFile(outNameFile);
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "DeconvoluteCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
