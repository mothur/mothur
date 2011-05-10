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

//**********************************************************************************************************************
vector<string> DegapSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DegapSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The degap.seqs command reads a fastafile and removes all gap characters.\n";
		helpString += "The degap.seqs command parameter is fasta.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your sequences, and is required unless you have a valid current fasta file. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The degap.seqs command should be in the following format: \n";
		helpString += "degap.seqs(fasta=yourFastaFile) \n";	
		helpString += "Example: degap.seqs(fasta=abrecovery.align) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
DegapSeqsCommand::DegapSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "DegapSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************
DegapSeqsCommand::DegapSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", false);
			if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { fastaFileNames.push_back(fastafile); m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else { 
				m->splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = m->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
						if (inputDir != "") {
							string path = m->hasPath(fastaFileNames[i]);
							//if the user has not given a path then, add inputdir. else leave path alone.
							if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
						}
		
						ifstream in;
						int ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
					
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getDefaultPath() != "") { //default path is set
								string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
								m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
								fastaFileNames[i] = tryPath;
							}
						}
						
						//if you can't open it, try default location
						if (ableToOpen == 1) {
							if (m->getOutputDir() != "") { //default path is set
								string tryPath = m->getOutputDir() + m->getSimpleName(fastaFileNames[i]);
								m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
								ifstream in2;
								ableToOpen = m->openInputFile(tryPath, in2, "noerror");
								in2.close();
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
				}
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

		}
	}
	catch(exception& e) {
		m->errorOut(e, "DegapSeqsCommand", "DegapSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************
int DegapSeqsCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Degapping sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
			ifstream inFASTA;
			m->openInputFile(fastaFileNames[s], inFASTA);
			
			ofstream outFASTA;
			string tempOutputDir = outputDir;
			if (outputDir == "") { tempOutputDir = m->hasPath(fastaFileNames[s]); }
			string degapFile = tempOutputDir + m->getRootName(m->getSimpleName(fastaFileNames[s])) + "ng.fasta";
			m->openOutputFile(degapFile, outFASTA);
			
			while(!inFASTA.eof()){
				if (m->control_pressed) {   outputTypes.clear(); inFASTA.close();  outFASTA.close(); remove(degapFile.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} return 0; }
				 
				Sequence currSeq(inFASTA);  m->gobble(inFASTA);
				if (currSeq.getName() != "") {
					outFASTA << ">" << currSeq.getName() << endl;
					outFASTA << currSeq.getUnaligned() << endl;
				}
			}
			inFASTA.close();
			outFASTA.close();
			
			outputNames.push_back(degapFile); outputTypes["fasta"].push_back(degapFile);
			
			if (m->control_pressed) {  outputTypes.clear(); remove(degapFile.c_str()); for (int j = 0; j < outputNames.size(); j++) {	remove(outputNames[j].c_str());	} return 0; }
		}
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
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

