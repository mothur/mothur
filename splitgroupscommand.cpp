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
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n";
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
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
		
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {			
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") {  groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { 			
				groupfile = m->getGroupFile(); 
				if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current groupfile and the group parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);		
			if (groups == "not found") { groups = ""; }
			else { m->splitAtDash(groups, Groups);	}
						
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(groupfile);	}
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
		
		groupMap = new GroupMap(groupfile);
		groupMap->readMap();
		
		SharedUtil util;  util.setGroups(Groups, groupMap->namesOfGroups);  
		
		if (namefile != "") {  readNames();  }
		splitFasta();
		
		delete groupMap;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} return 0; }
		
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
/**********************************************************************************************************************/
int SplitGroupCommand::readNames() { 
	try {
		//open input file
		ifstream in;
		m->openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol >> secondCol; m->gobble(in);
			
			vector<string> names;
			m->splitAtComma(secondCol, names);
			
			nameMap[firstCol] = names;
		}
		in.close();
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "readNames");
		exit(1);
	}
}

/**********************************************************************************************************************/
int SplitGroupCommand::splitFasta() { 
	try {
		
		string filerootFasta =  outputDir + m->getRootName(m->getSimpleName(fastafile));
		string filerootName =  outputDir + m->getRootName(m->getSimpleName(namefile));
		ofstream* temp;
		ofstream* temp2;
		map<string, ofstream*> filehandles;
		map<string, ofstream*>::iterator it3;

		for (int i=0; i<Groups.size(); i++) {
			temp = new ofstream;
			filehandles[Groups[i]+"fasta"] = temp;
			m->openOutputFile(filerootFasta + Groups[i] + ".fasta", *(filehandles[Groups[i]+"fasta"]));
			outputNames.push_back(filerootFasta + Groups[i] + ".fasta"); outputTypes["fasta"].push_back(filerootFasta + Groups[i] + ".fasta");
			
			if (namefile != "") {
				temp2 = new ofstream;
				filehandles[Groups[i]+"name"] = temp2;
				m->openOutputFile(filerootName + Groups[i] + ".names", *(filehandles[Groups[i]+"name"]));
				outputNames.push_back(filerootName + Groups[i] + ".names"); outputTypes["name"].push_back(filerootFasta + Groups[i] + ".names");
			}
		}
			
		//open input file
		ifstream in;
		m->openInputFile(fastafile, in);
	
		while (!in.eof()) {
			if (m->control_pressed) { break; }
		
			Sequence seq(in); m->gobble(in);
				
			if (seq.getName() != "") { 
				
				itNames = nameMap.find(seq.getName());
				
				if (itNames == nameMap.end()) {
					if (namefile != "") {  m->mothurOut(seq.getName() + " is not in your namesfile, ignoring."); m->mothurOutEndLine();  }
					else { //no names file given
						string group = groupMap->getGroup(seq.getName());
						
						if (m->inUsersGroups(group, Groups)) { //only add if this is in a group we want
							seq.printSequence(*(filehandles[group+"fasta"]));
						}else if(group == "not found") {
							m->mothurOut(seq.getName() + " is not in your groupfile. Ignoring."); m->mothurOutEndLine();
						}
					}
				}else{
					vector<string> names = itNames->second;
					map<string, string> group2Names;
					map<string, string>::iterator it;
					
					for (int i = 0; i < names.size(); i++) {  //build strings for new group names file, will select rep below
						string group = groupMap->getGroup(names[i]);
						
						if (m->inUsersGroups(group, Groups)) { //only add if this is in a group we want
							it = group2Names.find(group);
							if (it == group2Names.end()) {
								group2Names[group] = names[i] + ",";
							}else{
								group2Names[group] += names[i] + ",";
							}
						}else if(group == "not found") {
							m->mothurOut(names[i] + " is not in your groupfile. Ignoring."); m->mothurOutEndLine();
						}
					}
				
					for (map<string, string>::iterator itGroups = group2Names.begin(); itGroups != group2Names.end(); itGroups++) {
						//edit names string, by grabbing the first guy as the rep and removing the last comma
						string newNames = itGroups->second;
						newNames = newNames.substr(0, newNames.length()-1); 
						string repName = "";
						
						int pos = newNames.find_first_of(',');
						if (pos == newNames.npos) { //only one sequence so it represents itself
							repName = newNames;
						}else{
							repName = newNames.substr(0, pos);
						}
						
						*(filehandles[itGroups->first+"name"]) << repName << '\t' << newNames << endl;
						*(filehandles[itGroups->first+"fasta"]) << ">" << repName << endl << seq.getAligned() << endl;
					}
				}
			}
		}
			
		in.close();
			
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { 
			(*(filehandles[it3->first])).close();
			delete it3->second;
		}
		
		vector<string> newOutputNames;
		//remove blank files
		for (int i = 0; i < outputNames.size(); i++) {
			if (m->isBlank(outputNames[i])) {
				remove(outputNames[i].c_str());
			}else { newOutputNames.push_back(outputNames[i]); }
		}
		outputNames = newOutputNames;
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "splitFasta");
		exit(1);
	}
}
/**********************************************************************************************************************/

