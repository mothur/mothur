/*
 *  countgroupscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 8/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "countgroupscommand.h"
#include "sharedutilities.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> CountGroupsCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "sharedGroup", "sharedGroup", "none",false,false); parameters.push_back(pshared);
		CommandParameter pgroup("group", "InputTypes", "", "", "sharedGroup", "sharedGroup", "none",false,false); parameters.push_back(pgroup);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(paccnos);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountGroupsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The count.groups command counts sequences from a specfic group or set of groups from the following file types: group or shared file.\n";
		helpString += "The count.groups command parameters are accnos, group, shared and groups. You must provide a group or shared file.\n";
		helpString += "The accnos parameter allows you to provide a file containing the list of groups.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like.  You can separate group names with dashes.\n";
		helpString += "The count.groups command should be in the following format: count.groups(accnos=yourAccnos, group=yourGroupFile).\n";
		helpString += "Example count.groups(accnos=amazon.accnos, group=amazon.groups).\n";
		helpString += "or count.groups(groups=pasture, group=amazon.groups).\n";
		helpString += "Note: No spaces between parameter labels (i.e. group), '=' and parameters (i.e.yourGroupFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
CountGroupsCommand::CountGroupsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "CountGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
CountGroupsCommand::CountGroupsCommand(string option)  {
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
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  accnosfile = ""; }
			else { m->setAccnosFile(accnosfile); }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else {
				m->splitAtDash(groups, Groups);
				m->Groups = Groups;
			}
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
			else { m->setSharedFile(sharedfile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") {  	groupfile = "";	}
			else { m->setGroupFile(groupfile); }	
			
			if ((sharedfile == "") && (groupfile == "")) { 
				//give priority to shared, then group
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					groupfile = m->getGroupFile(); 
					if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You have no current groupfile or sharedfile and one is required."); m->mothurOutEndLine(); abort = true;
					}
				}
			}
			
			if ((accnosfile == "") && (Groups.size() == 0)) { Groups.push_back("all"); m->Groups = Groups; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "CountGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountGroupsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//get groups you want to remove
		if (accnosfile != "") { readAccnos(); }
		
		if (groupfile != "") {
			GroupMap groupMap(groupfile);
			groupMap.readMap();
			
			//make sure groups are valid
			//takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil util;
			util.setGroups(Groups, groupMap.namesOfGroups);
			
			for (int i = 0; i < Groups.size(); i++) {
				m->mothurOut(Groups[i] + " contains " + toString(groupMap.getNumSeqs(Groups[i])) + "."); m->mothurOutEndLine();
			}
		}
		
		if (m->control_pressed) { return 0; }
		
		if (sharedfile != "")		{		
			InputData input(sharedfile, "sharedfile");
			vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
			
			for (int i = 0; i < lookup.size(); i++) {
				m->mothurOut(lookup[i]->getGroup() + " contains " + toString(lookup[i]->getNumSeqs()) + "."); m->mothurOutEndLine();
				delete lookup[i];
			}			
		}
				
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void CountGroupsCommand::readAccnos(){
	try {
		Groups.clear();
		
		ifstream in;
		m->openInputFile(accnosfile, in);
		string name;
		
		while(!in.eof()){
			in >> name;
			
			Groups.push_back(name);
			
			m->gobble(in);
		}
		in.close();		
		
		m->Groups = Groups;
		
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "readAccnos");
		exit(1);
	}
}
//**********************************************************************************************************************


