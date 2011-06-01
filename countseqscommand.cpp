/*
 *  countseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/1/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "countseqscommand.h"
#include "groupmap.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> CountSeqsCommand::setParameters(){	
	try {
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The count.seqs command reads a name file and outputs a .seq.count file.  You may also provide a group file to get the counts broken down by group.\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include in the counts, by default all groups in your groupfile are used.\n";
		helpString += "When you use the groups parameter and a sequence does not represent any sequences from the groups you specify it is not included in the .count.summary file.\n";
		helpString += "The count.seqs command should be in the following format: count.seqs(name=yourNameFile).\n";
		helpString += "Example count.seqs(name=amazon.names).\n";
		helpString += "Note: No spaces between parameter labels (i.e. name), '=' and parameters (i.e.yourNameFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
CountSeqsCommand::CountSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "CountSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

CountSeqsCommand::CountSeqsCommand(string option)  {
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
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){					
				namefile = m->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current namefile and the name parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "all"; }
			m->splitAtDash(groups, Groups);
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(namefile);		}

		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "CountSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		ofstream out;
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(namefile)) + ".seq.count";
		m->openOutputFile(outputFileName, out); outputTypes["summary"].push_back(outputFileName);
		out << "Representative Sequence\t total\t";
		
		GroupMap* groupMap;
		if (groupfile != "") { 
			groupMap = new GroupMap(groupfile); groupMap->readMap(); 
			
			//make sure groups are valid. takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			util->setGroups(Groups, groupMap->namesOfGroups);
			delete util;
			
			//sort groupNames so that the group title match the counts below, this is needed because the map object automatically sorts
			sort(Groups.begin(), Groups.end());
			
			//print groupNames
			for (int i = 0; i < Groups.size(); i++) {
				out << Groups[i] << '\t';
			}
		}
		out << endl;
		
		//open input file
		ifstream in;
		m->openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol >> secondCol; m->gobble(in);
			
			vector<string> names;
			m->splitAtChar(secondCol, names, ',');
			
			if (groupfile != "") {
				//set to 0
				map<string, int> groupCounts;
				int total = 0;
				for (int i = 0; i < Groups.size(); i++) { groupCounts[Groups[i]] = 0; }
				
				//get counts for each of the users groups
				for (int i = 0; i < names.size(); i++) {
					string group = groupMap->getGroup(names[i]);
					
					if (group == "not found") { m->mothurOut("[ERROR]: " + names[i] + " is not in your groupfile, please correct."); m->mothurOutEndLine(); }
					else {
						map<string, int>::iterator it = groupCounts.find(group);
						
						//if not found, then this sequence is not from a group we care about
						if (it != groupCounts.end()) {
							it->second++;
							total++;
						}
					}
				}
				
				if (total != 0) {
					out << firstCol << '\t' << total << '\t';
					for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
						out << it->second << '\t';
					}
					out << endl;
				}
			}else {
				out << firstCol << '\t' << names.size() << endl;
			}
			
			
		}
		in.close();
		
		if (groupfile != "") { delete groupMap; }
		
		if (m->control_pressed) { remove(outputFileName.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "CountSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
