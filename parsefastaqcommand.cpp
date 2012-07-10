/*
 *  parsefastaqcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "parsefastaqcommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> ParseFastaQCommand::setParameters(){	
	try {
		CommandParameter pfastq("fastq", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfastq);
		CommandParameter pfasta("fasta", "Bool", "", "T", "", "", "",false,false); parameters.push_back(pfasta);
		CommandParameter pqual("qfile", "Bool", "", "T", "", "", "",false,false); parameters.push_back(pqual);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParseFastaQCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The fastq.info command reads a fastq file and creates a fasta and quality file.\n";
		helpString += "The fastq.info command parameters are fastq, fasta and qfile; fastq is required.\n";
		helpString += "The fastq.info command should be in the following format: fastq.info(fastaq=yourFastaQFile).\n";
        helpString += "The fasta parameter allows you to indicate whether you want a fasta file generated. Default=T.\n";
        helpString += "The qfile parameter allows you to indicate whether you want a quality file generated. Default=T.\n";
		helpString += "Example fastq.info(fastaq=test.fastaq).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fastq), '=' and yourFastQFile.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParseFastaQCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "fasta") {  outputFileName =  "fasta"; }
            else if (type == "qfile") {  outputFileName =  "qual"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getOutputFileNameTag");
		exit(1);
	}
}

//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(string option){
	try {
		abort = false; calledHelp = false;   
		
		if(option == "help") {	help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;

			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
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
				it = parameters.find("fastq");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastq"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			fastaQFile = validParameter.validFile(parameters, "fastq", true);
			if (fastaQFile == "not found") {	m->mothurOut("fastq is a required parameter for the fastq.info command.");	m->mothurOutEndLine();	abort = true;	}
			else if (fastaQFile == "not open")	{	fastaQFile = ""; abort = true;	}	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);	if (outputDir == "not found"){	outputDir = m->hasPath(fastaQFile); 	}
			
			string temp;
			temp = validParameter.validFile(parameters, "fasta", false);	if(temp == "not found"){	temp = "T";	}
			fasta = m->isTrue(temp); 

			temp = validParameter.validFile(parameters, "qfile", false);	if(temp == "not found"){	temp = "T";	}
			qual = m->isTrue(temp); 
			
			if ((!fasta) && (!qual)) { m->mothurOut("[ERROR]: no outputs selected. Aborting."); m->mothurOutEndLine(); abort=true; }

		}		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ParseFastaQCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//open Output Files
		string fastaFile = outputDir + m->getRootName(m->getSimpleName(fastaQFile)) + getOutputFileNameTag("fasta");
		string qualFile = outputDir + m->getRootName(m->getSimpleName(fastaQFile)) + getOutputFileNameTag("qfile");
		ofstream outFasta, outQual;
		
		if (fasta) { m->openOutputFile(fastaFile, outFasta);  outputNames.push_back(fastaFile); outputTypes["fasta"].push_back(fastaFile);	}
		if (qual) { m->openOutputFile(qualFile, outQual);	outputNames.push_back(qualFile);  outputTypes["qfile"].push_back(qualFile);		}
		
		ifstream in;
		m->openInputFile(fastaQFile, in);
		
		while (!in.eof()) {
			
			if (m->control_pressed) { break; }
		
			//read sequence name
			string name = m->getline(in); m->gobble(in);
			if (name == "") {  m->mothurOut("[ERROR]: Blank fasta name."); m->mothurOutEndLine(); m->control_pressed = true; break; }
			else if (name[0] != '@') { m->mothurOut("[ERROR]: reading " + name + " expected a name with @ as a leading character."); m->mothurOutEndLine(); m->control_pressed = true; break; }
			else { name = name.substr(1); }
			
			//read sequence
			string sequence = m->getline(in); m->gobble(in);
			if (sequence == "") {  m->mothurOut("[ERROR]: missing sequence for " + name); m->mothurOutEndLine(); m->control_pressed = true; break; }
			
			//read sequence name
			string name2 = m->getline(in); m->gobble(in);
			if (name2 == "") {  m->mothurOut("[ERROR]: Blank quality name."); m->mothurOutEndLine(); m->control_pressed = true; break; }
			else if (name2[0] != '+') { m->mothurOut("[ERROR]: reading " + name2 + " expected a name with + as a leading character."); m->mothurOutEndLine(); m->control_pressed = true; break; }
			else { name2 = name2.substr(1);  }
			
			//read quality scores
			string quality = m->getline(in); m->gobble(in);
			if (quality == "") {  m->mothurOut("[ERROR]: missing quality for " + name2); m->mothurOutEndLine(); m->control_pressed = true; break; }
			
			//sanity check sequence length and number of quality scores match
			if (name2 != "") { if (name != name2) { m->mothurOut("[ERROR]: names do not match. read " + name + " for fasta and " + name2 + " for quality."); m->mothurOutEndLine(); m->control_pressed = true; break; } }
			if (quality.length() != sequence.length()) { m->mothurOut("[ERROR]: Lengths do not match for sequence " + name + ". Read " + toString(sequence.length()) + " characters for fasta and " + toString(quality.length()) + " characters for quality scores."); m->mothurOutEndLine(); m->control_pressed = true; break; }
			
			//print sequence info to files
			if (fasta) { outFasta << ">" << name << endl << sequence << endl; }
			
			if (qual) { 
				vector<int> qualScores = convertQual(quality);
				outQual << ">" << name << endl;
				for (int i = 0; i < qualScores.size(); i++) { outQual << qualScores[i] << " "; }
				outQual << endl;
			}
		}
		
		in.close();
		if (fasta)	{ outFasta.close();	}
		if (qual)	{ outQual.close();	}
		
		if (m->control_pressed) { outputTypes.clear(); m->mothurRemove(fastaFile); m->mothurRemove(qualFile); return 0; }
		
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
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<int> ParseFastaQCommand::convertQual(string qual) {
	try {
		vector<int> qualScores;
		
		int controlChar = int('@');
		
		for (int i = 0; i < qual.length(); i++) { 
			int temp = int(qual[i]);
			temp -= controlChar;
			
			qualScores.push_back(temp);
		}
		
		return qualScores;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "convertQual");
		exit(1);
	}
}
//**********************************************************************************************************************



