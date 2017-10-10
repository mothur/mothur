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
vector<string> ReverseSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "fastaQual", "none","fasta",false,false,true); parameters.push_back(pfasta);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "fastaQual", "none","qfile",false,false,true); parameters.push_back(pqfile);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ReverseSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The reverse.seqs command reads a fastafile and outputs a fasta file containing the reverse compliment.\n";
		helpString += "The reverse.seqs command parameters fasta or qfile are required.\n";
		helpString += "The reverse.seqs command should be in the following format: \n";
		helpString += "reverse.seqs(fasta=yourFastaFile) \n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ReverseSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],rc,[extension]"; } 
        else if (type == "qfile") {  pattern = "[filename],rc,[extension]"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ReverseSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ReverseSeqsCommand::ReverseSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "ReverseSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

ReverseSeqsCommand::ReverseSeqsCommand(string option)  {
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
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastaFileName = validParameter.validFile(parameters, "fasta", true);
			if (fastaFileName == "not open") { abort = true; }
			else if (fastaFileName == "not found") { fastaFileName = "";}// m->mothurOut("fasta is a required parameter for the reverse.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			else { m->setFastaFile(fastaFileName); }
			
			qualFileName = validParameter.validFile(parameters, "qfile", true);
			if (qualFileName == "not open") { abort = true; }
			else if (qualFileName == "not found") { qualFileName = ""; }//m->mothurOut("fasta is a required parameter for the reverse.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			else { m->setQualFile(qualFileName); }
			
			if((fastaFileName == "") && (qualFileName == "")){
				fastaFileName = m->getFastaFile(); 
				if (fastaFileName != "") {  m->mothurOut("Using " + fastaFileName + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 
					qualFileName = m->getQualFile(); 
					if (qualFileName != "") {  m->mothurOut("Using " + qualFileName + " as input file for the qfile parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You have no current files for fasta or qfile, and fasta or qfile is a required parameter for the reverse.seqs command."); m->mothurOutEndLine();
						abort = true;
					}
				}
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
			}

		}
	}
	catch(exception& e) {
		m->errorOut(e, "ReverseSeqsCommand", "ReverseSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************


int ReverseSeqsCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		string fastaReverseFileName;
		
		if(fastaFileName != ""){
			ifstream inFASTA;
			m->openInputFile(fastaFileName, inFASTA);
			
			ofstream outFASTA;
			string tempOutputDir = outputDir;
			if (outputDir == "") { tempOutputDir += m->hasPath(fastaFileName); } //if user entered a file with a path then preserve it
            map<string, string> variables; 
            variables["[filename]"] = tempOutputDir + m->getRootName(m->getSimpleName(fastaFileName));
            variables["[extension]"] = m->getExtension(fastaFileName);
			fastaReverseFileName = getOutputFileName("fasta", variables);
			m->openOutputFile(fastaReverseFileName, outFASTA);
			
			while(!inFASTA.eof()){
				if (m->getControl_pressed()) {  inFASTA.close();  outFASTA.close(); m->mothurRemove(fastaReverseFileName); return 0; }
				 
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
			string tempOutputDir = outputDir;
			if (outputDir == "") { tempOutputDir += m->hasPath(qualFileName); } //if user entered a file with a path then preserve it
            map<string, string> variables; 
            variables["[filename]"] = tempOutputDir + m->getRootName(m->getSimpleName(qualFileName));
            variables["[extension]"] = m->getExtension(qualFileName);
			string qualReverseFileName = getOutputFileName("qfile", variables);
            m->openOutputFile(qualReverseFileName, outQual);

			while(!inQual.eof()){
				if (m->getControl_pressed()) {  inQual.close();  outQual.close(); m->mothurRemove(qualReverseFileName); return 0; }
				currQual = QualityScores(inQual);  m->gobble(inQual);
				currQual.flipQScores();	
				currQual.printQScores(outQual);
			}
			inQual.close();
			outQual.close();
			outputNames.push_back(qualReverseFileName); outputTypes["qfile"].push_back(qualReverseFileName);
		}
		
		if (m->getControl_pressed()) {  m->mothurRemove(qualReverseFileName); m->mothurRemove(fastaReverseFileName); return 0; }
		
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
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
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
