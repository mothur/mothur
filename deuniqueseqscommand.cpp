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
vector<string> DeUniqueSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pname);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DeUniqueSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The deunique.seqs command reads a fastafile and namefile, and creates a fastafile containing all the sequences.\n";
		helpString += "The deunique.seqs command parameters are fasta and name, both are required, unless you have valid current name and fasta files.\n";
		helpString += "The deunique.seqs command should be in the following format: \n";
		helpString += "deunique.seqs(fasta=yourFastaFile, name=yourNameFile) \n";	
		helpString += "Example deunique.seqs(fasta=abrecovery.unique.fasta, name=abrecovery.names).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
DeUniqueSeqsCommand::DeUniqueSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "DeconvoluteCommand");
		exit(1);
	}
}
/**************************************************************************************/
DeUniqueSeqsCommand::DeUniqueSeqsCommand(string option)  {	
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			vector<string> myArray = setParameters();
			
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
			else if (fastaFile == "not found") { 				
				fastaFile = m->getFastaFile(); 
				if (fastaFile != "") { m->mothurOut("Using " + fastaFile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastaFile); //if user entered a file with a path then preserve it	
			}
			
			nameFile = validParameter.validFile(parameters, "name", true);
			if (nameFile == "not open") { abort = true; }
			else if (nameFile == "not found"){					
				nameFile = m->getNameFile(); 
				if (nameFile != "") { m->mothurOut("Using " + nameFile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current namefile and the name parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "DeUniqueSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
int DeUniqueSeqsCommand::execute() {	
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

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
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
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
