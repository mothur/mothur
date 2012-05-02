/*
 *  splitgroupscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitgroupscommand.h"
#include "sharedutilities.h"
#include "sequenceparser.h"

//**********************************************************************************************************************
vector<string> SplitGroupCommand::setParameters(){	
	try {		
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SplitGroupCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The split.group command reads a group file, and parses your fasta and names files by groups. \n";
		helpString += "The split.group command parameters are fasta, name, group and groups.\n";
		helpString += "The fasta and group parameters are required.\n";
		helpString += "The groups parameter allows you to select groups to create files for.  \n";
		helpString += "For example if you set groups=A-B-C, you will get a .A.fasta, .A.names, .B.fasta, .B.names, .C.fasta, .C.names files.  \n";
		helpString += "If you want .fasta and .names files for all groups, set groups=all.  \n";
		helpString += "The split.group command should be used in the following format: split.group(fasta=yourFasta, group=yourGroupFile).\n";
		helpString += "Example: split.group(fasta=abrecovery.fasta, group=abrecovery.groups).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
SplitGroupCommand::SplitGroupCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "SplitGroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SplitGroupCommand::SplitGroupCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
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

			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			else { m->setNameFile(namefile); }
		
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {			
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); }	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") {  groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { 			
				groupfile = m->getGroupFile(); 
				if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current groupfile and the group parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else {  m->setGroupFile(groupfile); }
			
			groups = validParameter.validFile(parameters, "groups", false);		
			if (groups == "not found") { groups = ""; }
			else { m->splitAtDash(groups, Groups);	}
						
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(groupfile);	}
			
			if (namefile == "") {
				vector<string> files; files.push_back(fastafile);
				parser.getNameFile(files);
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "SplitAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SplitGroupCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		SequenceParser* parser;
		
		if (namefile == "") {	parser = new SequenceParser(groupfile, fastafile);				}
		else				{	parser = new SequenceParser(groupfile, fastafile, namefile);	}
		
		if (m->control_pressed) { delete parser; return 0; }

		vector<string> namesGroups = parser->getNamesOfGroups();
		SharedUtil util;  util.setGroups(Groups, namesGroups);  
		
		string fastafileRoot = outputDir + m->getRootName(m->getSimpleName(fastafile));
		string namefileRoot = outputDir + m->getRootName(m->getSimpleName(namefile));
		
		m->mothurOutEndLine();
		for (int i = 0; i < Groups.size(); i++) {
			
			m->mothurOut("Processing group: " + Groups[i]); m->mothurOutEndLine();
			
			string newFasta = fastafileRoot + Groups[i] + ".fasta";
			string newName = namefileRoot + Groups[i] + ".names";
			
			parser->getSeqs(Groups[i], newFasta, false);
			outputNames.push_back(newFasta); outputTypes["fasta"].push_back(newFasta);
			if (m->control_pressed) { delete parser; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} return 0; }

			if (namefile != "") { 
				parser->getNameMap(Groups[i], newName); 
				outputNames.push_back(newName); outputTypes["name"].push_back(newName);
			}
			
			if (m->control_pressed) { delete parser; for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} return 0; }
		}
		
		delete parser;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);	} return 0; }
		
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


