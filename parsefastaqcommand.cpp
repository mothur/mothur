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
vector<string> ParseFastaQCommand::getValidParameters(){	
	try {
		string Array[] =  {"fastq", "outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qual"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ParseFastaQCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fastq"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ParseFastaQCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
ParseFastaQCommand::ParseFastaQCommand(string option){
	try {
		abort = false;
		
		if(option == "help") {	help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fastq", "outputdir", "inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
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
			outputTypes["qual"] = tempOutNames;
			
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

		}		
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "ParseFastaQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ParseFastaQCommand::help(){
	try {
		m->mothurOut("The fastq.info command reads a fastq file and creates a fasta and quality file.\n");
		m->mothurOut("The fastq.info command parameter is fastq, and it is required.\n");
		m->mothurOut("The fastq.info command should be in the following format: fastq.info(fastaq=yourFastaQFile).\n");
		m->mothurOut("Example fastq.info(fastaq=test.fastaq).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fastq), '=' and yourFastQFile.\n");
		m->mothurOutEndLine();
	}
	catch(exception& e) {
		m->errorOut(e, "ParseFastaQCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************

ParseFastaQCommand::~ParseFastaQCommand()	{	/*	do nothing	*/	}

//**********************************************************************************************************************

int ParseFastaQCommand::execute(){
	try {
		if (abort == true) {	return 0;	}
		
		//open Output Files
		string fastaFile = outputDir + m->getRootName(m->getSimpleName(fastaQFile)) + "fasta";
		string qualFile = outputDir + m->getRootName(m->getSimpleName(fastaQFile)) + "qual";
		ofstream outFasta, outQual;
		m->openOutputFile(fastaFile, outFasta);  outputNames.push_back(fastaFile); outputTypes["fasta"].push_back(fastaFile);
		m->openOutputFile(qualFile, outQual);	outputNames.push_back(qualFile);  outputTypes["qual"].push_back(qualFile);
		
		ifstream in;
		m->openInputFile(fastaQFile, in);
		
		while (!in.eof()) {
		
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
			string qual = m->getline(in); m->gobble(in);
			if (qual == "") {  m->mothurOut("[ERROR]: missing quality for " + name2); m->mothurOutEndLine(); m->control_pressed = true; break; }
			
			//sanity check sequence length and number of quality scores match
			if (name2 != "") { if (name != name2) { m->mothurOut("[ERROR]: names do not match. read " + name + " for fasta and " + name2 + " for quality."); m->mothurOutEndLine(); m->control_pressed = true; break; } }
			if (qual.length() != sequence.length()) { m->mothurOut("[ERROR]: lengths do not match. read " + toString(sequence.length()) + " characters for fasta and " + toString(qual.length()) + " characters for quality scores."); m->mothurOutEndLine(); m->control_pressed = true; break; }
			
			//convert quality scores
			vector<int> qualScores = convertQual(qual);
			
			//print sequence info to files
			outFasta << ">" << name << endl << sequence << endl;
			
			outQual << ">" << name << endl;
			for (int i = 0; i < qualScores.size(); i++) { outQual << qualScores[i] << " "; }
			outQual << endl;
		}
		
		in.close();
		outFasta.close();
		outQual.close();
		
		if (m->control_pressed) { outputTypes.clear(); remove(fastaFile.c_str()); remove(qualFile.c_str()); return 0; }
		
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
		
		int controlChar = int('!');
		
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



