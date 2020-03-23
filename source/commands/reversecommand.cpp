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
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["qfile"] = tempOutNames;
		
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
//***************************************************************************************************************

ReverseSeqsCommand::ReverseSeqsCommand(string option)  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastaFileName = validParameter.validFile(parameters, "fasta");
			if (fastaFileName == "not open") { abort = true; }
			else if (fastaFileName == "not found") { fastaFileName = "";}// m->mothurOut("fasta is a required parameter for the reverse.seqs command.\n");  abort = true;  }	
			else { current->setFastaFile(fastaFileName); }
			
			qualFileName = validParameter.validFile(parameters, "qfile");
			if (qualFileName == "not open") { abort = true; }
			else if (qualFileName == "not found") { qualFileName = ""; }//m->mothurOut("fasta is a required parameter for the reverse.seqs command.\n");  abort = true;  }	
			else { current->setQualFile(qualFileName); }
			
			if((fastaFileName == "") && (qualFileName == "")){
				fastaFileName = current->getFastaFile(); 
				if (fastaFileName != "") {  m->mothurOut("Using " + fastaFileName + " as input file for the fasta parameter.\n");  }
				else { 
					qualFileName = current->getQualFile(); 
					if (qualFileName != "") {  m->mothurOut("Using " + qualFileName + " as input file for the qfile parameter.\n");  }
					else { 
						m->mothurOut("You have no current files for fasta or qfile, and fasta or qfile is a required parameter for the reverse.seqs command.\n"); 
						abort = true;
					}
				}
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
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
			util.openInputFile(fastaFileName, inFASTA);
			
			ofstream outFASTA;
			string tempOutputDir = outputDir;
			if (outputDir == "") { tempOutputDir += util.hasPath(fastaFileName); } //if user entered a file with a path then preserve it
            map<string, string> variables; 
            variables["[filename]"] = tempOutputDir + util.getRootName(util.getSimpleName(fastaFileName));
            variables["[extension]"] = util.getExtension(fastaFileName);
			fastaReverseFileName = getOutputFileName("fasta", variables);
			util.openOutputFile(fastaReverseFileName, outFASTA);
			
			while(!inFASTA.eof()){
				if (m->getControl_pressed()) {  inFASTA.close();  outFASTA.close(); util.mothurRemove(fastaReverseFileName); return 0; }
				 
				Sequence currSeq(inFASTA);  util.gobble(inFASTA);
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
			util.openInputFile(qualFileName, inQual);
			
			ofstream outQual;
			string tempOutputDir = outputDir;
			if (outputDir == "") { tempOutputDir += util.hasPath(qualFileName); } //if user entered a file with a path then preserve it
            map<string, string> variables; 
            variables["[filename]"] = tempOutputDir + util.getRootName(util.getSimpleName(qualFileName));
            variables["[extension]"] = util.getExtension(qualFileName);
			string qualReverseFileName = getOutputFileName("qfile", variables);
            util.openOutputFile(qualReverseFileName, outQual);

			while(!inQual.eof()){
				if (m->getControl_pressed()) {  inQual.close();  outQual.close(); util.mothurRemove(qualReverseFileName); return 0; }
				currQual = QualityScores(inQual);  util.gobble(inQual);
				currQual.flipQScores();	
				currQual.printQScores(outQual);
			}
			inQual.close();
			outQual.close();
			outputNames.push_back(qualReverseFileName); outputTypes["qfile"].push_back(qualReverseFileName);
		}
		
		if (m->getControl_pressed()) {  util.mothurRemove(qualReverseFileName); util.mothurRemove(fastaReverseFileName); return 0; }
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
		}
		
		
		m->mothurOut("\nOutput File Names: \n"); 
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
