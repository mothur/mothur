/*
 *  degapseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "degapseqscommand.h"
#include "sequence.hpp"

//***************************************************************************************************************

DegapSeqsCommand::DegapSeqsCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the degap.seqs command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					if (inputDir != "") {
						string path = hasPath(fastaFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
					}
	
					ifstream in;
					int ableToOpen = openInputFile(fastaFileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ableToOpen = openInputFile(tryPath, in, "noerror");
							fastaFileNames[i] = tryPath;
						}
					}
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + fastaFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine();
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

		}
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "DegapSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void DegapSeqsCommand::help(){
	try {
		m->mothurOut("The degap.seqs command reads a fastafile and removes all gap characters.\n");
		m->mothurOut("The degap.seqs command parameter is fasta.\n");
		m->mothurOut("The fasta parameter allows you to enter the fasta file containing your sequences, and is required. \n");
		m->mothurOut("You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n");
		m->mothurOut("The degap.seqs command should be in the following format: \n");
		m->mothurOut("degap.seqs(fasta=yourFastaFile) \n");	
		m->mothurOut("Example: degap.seqs(fasta=abrecovery.align) \n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");	
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

DegapSeqsCommand::~DegapSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************


int DegapSeqsCommand::execute(){
	try{
		
		if (abort == true) { return 0; }
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Degapping sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			ifstream inFASTA;
			openInputFile(fastaFileNames[s], inFASTA);
			
			ofstream outFASTA;
			string degapFile = outputDir + getRootName(getSimpleName(fastaFileNames[s])) + "ng.fasta";
			openOutputFile(degapFile, outFASTA);
			
			while(!inFASTA.eof()){
				if (m->control_pressed) {  inFASTA.close();  outFASTA.close(); remove(degapFile.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} return 0; }
				 
				Sequence currSeq(inFASTA);  gobble(inFASTA);
				if (currSeq.getName() != "") {
					outFASTA << ">" << currSeq.getName() << endl;
					outFASTA << currSeq.getUnaligned() << endl;
				}
			}
			inFASTA.close();
			outFASTA.close();
			
			outputNames.push_back(degapFile);
			
			if (m->control_pressed) {  remove(degapFile.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} return 0; }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();

		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************

