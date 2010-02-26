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
DeconvoluteCommand::DeconvoluteCommand(string option)  {	
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "name","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			inFastaName = validParameter.validFile(parameters, "fasta", true);
			if (inFastaName == "not open") { abort = true; }
			else if (inFastaName == "not found") { inFastaName = ""; m->mothurOut("fasta is a required parameter for the unique.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(inFastaName); //if user entered a file with a path then preserve it	
			}
			
			oldNameMapFName = validParameter.validFile(parameters, "name", true);
			if (oldNameMapFName == "not open") { abort = true; }
			else if (oldNameMapFName == "not found"){	oldNameMapFName = "";	}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "DeconvoluteCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void DeconvoluteCommand::help(){
	try {
		m->mothurOut("The unique.seqs command reads a fastafile and creates a namesfile.\n");
		m->mothurOut("It creates a file where the first column is the groupname and the second column is a list of sequence names who have the same sequence. \n");
		m->mothurOut("If the sequence is unique the second column will just contain its name. \n");
		m->mothurOut("The unique.seqs command parameter is fasta and it is required.\n");
		m->mothurOut("The unique.seqs command should be in the following format: \n");
		m->mothurOut("unique.seqs(fasta=yourFastaFile) \n");	
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "help");
		exit(1);
	}
}

/**************************************************************************************/
int DeconvoluteCommand::execute() {	
	try {
		
		if (abort == true) { return 0; }

		//prepare filenames and open files
		string outNameFile = outputDir + getRootName(getSimpleName(inFastaName)) + "names";
		string outFastaFile = outputDir + getRootName(getSimpleName(inFastaName)) + "unique" + getExtension(inFastaName);
		
		FastaMap fastamap;
	
		if(oldNameMapFName == "")	{	fastamap.readFastaFile(inFastaName);					}
		else						{	fastamap.readFastaFile(inFastaName, oldNameMapFName);	}
		
		fastamap.printCondensedFasta(outFastaFile);
		fastamap.printNamesFile(outNameFile);
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outFastaFile); m->mothurOutEndLine();	
		m->mothurOut(outNameFile); m->mothurOutEndLine();
		m->mothurOutEndLine();


		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
