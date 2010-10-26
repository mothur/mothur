/*
 *  deuniqueseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "deuniqueseqscommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> DeUniqueSeqsCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta", "name","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
DeUniqueSeqsCommand::DeUniqueSeqsCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "DeconvoluteCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> DeUniqueSeqsCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta","name"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> DeUniqueSeqsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "getRequiredFiles");
		exit(1);
	}
}
/**************************************************************************************/
DeUniqueSeqsCommand::DeUniqueSeqsCommand(string option)  {	
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
		
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
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta", true);
			if (fastaFile == "not open") { abort = true; }
			else if (fastaFile == "not found") { fastaFile = ""; m->mothurOut("fasta is a required parameter for the deunique.seqs command."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastaFile); //if user entered a file with a path then preserve it	
			}
			
			nameFile = validParameter.validFile(parameters, "name", true);
			if (nameFile == "not open") { abort = true; }
			else if (nameFile == "not found"){	nameFile = "";	m->mothurOut("name is a required parameter for the deunique.seqs command."); m->mothurOutEndLine(); abort = true;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "DeUniqueSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void DeUniqueSeqsCommand::help(){
	try {
		m->mothurOut("The deunique.seqs command reads a fastafile and namefile, and creates a fastafile containing all the sequences.\n");
		m->mothurOut("The deunique.seqs command parameters are fasta and name, both are required.\n");
		m->mothurOut("The deunique.seqs command should be in the following format: \n");
		m->mothurOut("deunique.seqs(fasta=yourFastaFile, name=yourNameFile) \n");	
		m->mothurOut("Example deunique.seqs(fasta=abrecovery.unique.fasta, name=abrecovery.names).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "help");
		exit(1);
	}
}

/**************************************************************************************/
int DeUniqueSeqsCommand::execute() {	
	try {
		
		if (abort == true) { return 0; }

		//prepare filenames and open files
		ofstream out;
		string outFastaFile = m->getRootName(m->getSimpleName(fastaFile));
		int pos = outFastaFile.find("unique");
		if (pos != string::npos) {
			outFastaFile = outputDir + outFastaFile.substr(0, pos) + "redundant" + m->getExtension(fastaFile);
		}else{
			outFastaFile = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + "redundant" + m->getExtension(fastaFile);
		}
		m->openOutputFile(outFastaFile, out);
		
		readNamesFile();
		if (m->control_pressed) {  out.close(); outputTypes.clear(); remove(outFastaFile.c_str()); return 0; }
		
		ifstream in;
		m->openInputFile(fastaFile, in);
		
		while (!in.eof()) {
		
			if (m->control_pressed) { in.close(); out.close(); outputTypes.clear(); remove(outFastaFile.c_str()); return 0; }
			
			Sequence seq(in); m->gobble(in);
			
			if (seq.getName() != "") {
				
				//look for sequence name in nameMap
				map<string, string>::iterator it = nameMap.find(seq.getName());
				
				if (it == nameMap.end()) {	m->mothurOut("[ERROR]: Your namefile does not contain " + seq.getName() + ", aborting."); m->mothurOutEndLine(); m->control_pressed = true; }
				else {
					vector<string> names;
					m->splitAtComma(it->second, names);
					
					//output sequences
					for (int i = 0; i < names.size(); i++) {
						out << ">" << names[i] << endl;
						out << seq.getAligned() << endl;
					}
					
					//remove seq from name map so we can check for seqs in namefile not in fastafile later
					nameMap.erase(it);
				}
			}
		}
		in.close();
		out.close(); 
		
		if (nameMap.size() != 0) { //then there are names in the namefile not in the fastafile
			for (map<string, string>::iterator it = nameMap.begin(); it != nameMap.end(); it++) {  
				m->mothurOut(it->first + " is not in your fasta file, but is in your name file. Please correct."); m->mothurOutEndLine();
			}
		}
				
		if (m->control_pressed) { outputTypes.clear(); remove(outFastaFile.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outFastaFile); m->mothurOutEndLine();	
		outputNames.push_back(outFastaFile);  outputTypes["fasta"].push_back(outFastaFile);  
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int DeUniqueSeqsCommand::readNamesFile() {
	try {
		
		ifstream inNames;
		m->openInputFile(nameFile, inNames);
		
		string name, names;
		map<string, string>::iterator it;
	
		while(inNames){
			
			if(m->control_pressed) { break; }
			
			inNames >> name;	m->gobble(inNames);		
			inNames >> names;		
			
			it = nameMap.find(name);
			
			if (it == nameMap.end()) {	nameMap[name] = names; }
			else { m->mothurOut("[ERROR]: Your namefile already contains " + name + ", aborting."); m->mothurOutEndLine(); m->control_pressed = true; }
					
			m->gobble(inNames);
		}
		inNames.close();
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "readNamesFile");
		exit(1);
	}
}

/**************************************************************************************/
