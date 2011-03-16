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
#include "qualityscores.h"

//**********************************************************************************************************************
vector<string> ReverseSeqsCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta", "qfile", "outputdir", "inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ReverseSeqsCommand::ReverseSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "ReverseSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

vector<string> ReverseSeqsCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta", "qfile"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "getRequiredParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> ReverseSeqsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "getRequiredFiles");
		exit(1);
	}
}//***************************************************************************************************************

ReverseSeqsCommand::ReverseSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta", "qfile", "outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
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
			outputTypes["qfile"] = tempOutNames;

			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastaFileName = validParameter.validFile(parameters, "fasta", true);
			if (fastaFileName == "not open") { abort = true; }
			else if (fastaFileName == "not found") { fastaFileName = "";}// m->mothurOut("fasta is a required parameter for the reverse.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			
			qualFileName = validParameter.validFile(parameters, "qfile", true);
			if (qualFileName == "not open") { abort = true; }
			else if (qualFileName == "not found") { qualFileName = ""; }//m->mothurOut("fasta is a required parameter for the reverse.seqs command."); m->mothurOutEndLine(); abort = true;  }	

			if(fastaFileName == "not open" && qualFileName == "not open"){
				m->mothurOut("fasta or qfile is a required parameter for the reverse.seqs command.");
				m->mothurOutEndLine();
				abort = true;
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastaFileName); //if user entered a file with a path then preserve it	
			}

		}
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "ReverseSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ReverseSeqsCommand::help(){
	try {
		m->mothurOut("The reverse.seqs command reads a fastafile and outputs a fasta file containing the reverse compliment.\n");
		m->mothurOut("The reverse.seqs command parameters fasta or qfile are required.\n");
		m->mothurOut("The reverse.seqs command should be in the following format: \n");
		m->mothurOut("reverse.seqs(fasta=yourFastaFile) \n");	
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "help");
		exit(1);
	}
}

//***************************************************************************************************************

ReverseSeqsCommand::~ReverseSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************


int ReverseSeqsCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		string fastaReverseFileName;
		
		if(fastaFileName != ""){
			ifstream inFASTA;
			m->openInputFile(fastaFileName, inFASTA);
			
			ofstream outFASTA;
			fastaReverseFileName = outputDir + m->getRootName(m->getSimpleName(fastaFileName)) + "rc" + m->getExtension(fastaFileName);
			m->openOutputFile(fastaReverseFileName, outFASTA);
			
			while(!inFASTA.eof()){
				if (m->control_pressed) {  inFASTA.close();  outFASTA.close(); remove(fastaReverseFileName.c_str()); return 0; }
				 
				Sequence currSeq(inFASTA);  m->gobble(inFASTA);
				if (currSeq.getName() != "") {
					currSeq.reverseComplement();
					currSeq.printSequence(outFASTA);
				}
			}
			inFASTA.close();
			outFASTA.close();
			outputNames.push_back(fastaReverseFileName); outputTypes["fasta"].push_back(fastaReverseFileName);
		}
		
		string qualReverseFileName;

		if(qualFileName != ""){
			QualityScores currQual;

			ifstream inQual;
			m->openInputFile(qualFileName, inQual);
			
			ofstream outQual;
			string qualReverseFileName = outputDir + m->getRootName(m->getSimpleName(qualFileName)) + "rc" + m->getExtension(qualFileName);
			m->openOutputFile(qualReverseFileName, outQual);

			while(!inQual.eof()){
				if (m->control_pressed) {  inQual.close();  outQual.close(); remove(qualReverseFileName.c_str()); return 0; }
				currQual = QualityScores(inQual);  m->gobble(inQual);
				currQual.flipQScores();	
				currQual.printQScores(outQual);
			}
			inQual.close();
			outQual.close();
			outputNames.push_back(qualReverseFileName); outputTypes["qfile"].push_back(qualReverseFileName);
		}
		
		if (m->control_pressed) {  remove(qualReverseFileName.c_str()); remove(fastaReverseFileName.c_str()); return 0; }
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setQualFile(current); }
		}
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		for(int i=0;i<outputNames.size();i++){
			m->mothurOut(outputNames[i]);
			m->mothurOutEndLine();
		}
		
		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************
