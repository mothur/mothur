/*
 *  makefastqcommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "makefastqcommand.h"
#include "sequence.hpp"
#include "qualityscores.h"

//**********************************************************************************************************************
vector<string> MakeFastQCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta","qfile","outputdir","inputdir" };
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeFastQCommand::MakeFastQCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["fastq"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "MakeFastQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> MakeFastQCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta","qfile"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> MakeFastQCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeFastQCommand::MakeFastQCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","qfile", "outputdir","inputdir" };
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
			outputTypes["fastq"] = tempOutNames;
			
						
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
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; fastafile = ""; }
			else if (fastafile == "not found") {  fastafile = "";  m->mothurOut("You must provide a fasta file."); m->mothurOutEndLine(); abort = true; }	
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; qualfile = ""; }
			else if (qualfile == "not found") {  qualfile = "";  m->mothurOut("You must provide a quality file."); m->mothurOutEndLine(); abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(fastafile);		}

		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "MakeFastQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void MakeFastQCommand::help(){
	try {
		/*m->mothurOut("The get.groups command selects sequences from a specfic group or set of groups from the following file types: fasta, name, group, list, taxonomy.\n");
		m->mothurOut("It outputs a file containing the sequences in the those specified groups.\n");
		m->mothurOut("The get.groups command parameters are accnos, fasta, name, group, list, taxonomy and groups. The group parameter is required.\n");
		m->mothurOut("You must also provide an accnos containing the list of groups to get or set the groups parameter to the groups you wish to select.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like.  You can separate group names with dashes.\n");
		m->mothurOut("The get.groups command should be in the following format: get.groups(accnos=yourAccnos, fasta=yourFasta, group=yourGroupFile).\n");
		m->mothurOut("Example get.seqs(accnos=amazon.accnos, fasta=amazon.fasta, group=amazon.groups).\n");
		m->mothurOut("or get.seqs(groups=pasture, fasta=amazon.fasta, amazon.groups).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");*/
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int MakeFastQCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		
		string outputFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "fastq";
		outputNames.push_back(outputFile); outputTypes["fastq"].push_back(outputFile);
		
		ofstream out;
		m->openOutputFile(outputFile, out);
		
		ifstream qFile;
		m->openInputFile(qualfile, qFile);
		
		ifstream fFile;
		m->openInputFile(fastafile, fFile);
		
		while (!fFile.eof() && !qFile.eof()) {
			
			if (m->control_pressed) { break; }
			
			Sequence currSeq(fFile); m->gobble(fFile);
			QualityScores currQual(qFile);  m->gobble(qFile);
			
			if (currSeq.getName() != currQual.getName()) { m->mothurOut("[ERROR]: mismatch between fasta and quality files. Found " + currSeq.getName() + " in fasta file and " + currQual.getName() + " in quality file."); m->mothurOutEndLine(); m->control_pressed = true; }
			else {
				//print sequence
				out << '@' << currSeq.getName() << endl << currSeq.getAligned() << endl;
				
				string qualityString = convertQual(currQual.getQualityScores());
				
				//print quality info
				out << '+' << currQual.getName() << endl << qualityString << endl;
			}
			
		}
		
		fFile.close();
		qFile.close();
		out.close();
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
					
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeFastQCommand::convertQual(vector<int> qual) {
	try {
		string qualScores;
		
		int controlChar = int('!');
		
		for (int i = 0; i < qual.size(); i++) { 
			int temp = qual[i] + controlChar;
			char qualChar = (char) temp;
			
			qualScores += qualChar;
		}
		
		return qualScores;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "convertQual");
		exit(1);
	}
}
//**********************************************************************************************************************




		
