/*
 *  removegroupscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "removegroupscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"
#include "sharedutilities.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> RemoveGroupsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pfasta);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "FNGLT-sharedGroup", "none",false,false); parameters.push_back(pshared);
		CommandParameter pname("name", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "FNGLT-sharedGroup", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(ptaxonomy);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(paccnos);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveGroupsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.groups command removes sequences from a specfic group or set of groups from the following file types: fasta, name, group, list, taxonomy or sharedfile.\n";
		helpString += "It outputs a file containing the sequences NOT in the those specified groups, or with a sharedfile eliminates the groups you selected.\n";
		helpString += "The remove.groups command parameters are accnos, fasta, name, group, list, taxonomy, shared and groups. The group parameter is required, unless you have a current group file or are using a sharedfile.\n";
		helpString += "You must also provide an accnos containing the list of groups to remove or set the groups parameter to the groups you wish to remove.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like removed.  You can separate group names with dashes.\n";
		helpString += "The remove.groups command should be in the following format: remove.groups(accnos=yourAccnos, fasta=yourFasta, group=yourGroupFile).\n";
		helpString += "Example remove.groups(accnos=amazon.accnos, fasta=amazon.fasta, group=amazon.groups).\n";
		helpString += "or remove.groups(groups=pasture, fasta=amazon.fasta, amazon.groups).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
RemoveGroupsCommand::RemoveGroupsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "RemoveGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveGroupsCommand::RemoveGroupsCommand(string option)  {
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["shared"] = tempOutNames;
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
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
				
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
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
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
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
			if (accnosfile == "not open") { accnosfile = ""; abort = true; }
			else if (accnosfile == "not found") {  accnosfile = ""; }	
			else { m->setAccnosFile(accnosfile); }
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
			else { m->setFastaFile(fastafile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = "";  abort = true; }
			else if (groupfile == "not found") {  	groupfile = "";		}
			else { m->setGroupFile(groupfile); }	
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			else { m->setListFile(listfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not open") { taxfile = ""; abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			else { m->setTaxonomyFile(taxfile); }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
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
				//is there are current file available for any of these?
				if ((namefile != "") || (fastafile != "") || (listfile != "") || (taxfile != "")) {
					//give priority to group, then shared
					groupfile = m->getGroupFile(); 
					if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
					else { 
						sharedfile = m->getSharedFile(); 
						if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("You have no current groupfile or sharedfile and one is required."); m->mothurOutEndLine(); abort = true;
						}
					}
				}else {
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
			}
			
			if ((accnosfile == "") && (Groups.size() == 0)) { m->mothurOut("You must provide an accnos file containing group names or specify groups using the groups parameter."); m->mothurOutEndLine(); abort = true; }
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "")  && (sharedfile == "") && (listfile == "") && (taxfile == ""))  { m->mothurOut("You must provide at least one of the following: fasta, name, taxonomy, group, shared or list."); m->mothurOutEndLine(); abort = true; }
			if ((groupfile == "") && ((namefile != "") || (fastafile != "") || (listfile != "") || (taxfile != "")))  { m->mothurOut("If using a fasta, name, taxonomy, group or list, then you must provide a group file."); m->mothurOutEndLine(); abort = true; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "RemoveGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveGroupsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//get groups you want to remove
		if (accnosfile != "") { readAccnos(); }
		
		if (groupfile != "") {
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
			
			//make sure groups are valid
			//takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			vector<string> namesGroups = groupMap->getNamesOfGroups();
			util->setGroups(Groups, namesGroups);
			delete util;
			
			//fill names with names of sequences that are from the groups we want to remove 
			fillNames();
			
			delete groupMap;
		}
				
		if (m->control_pressed) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (namefile != "")			{		readName();		}
		if (fastafile != "")		{		readFasta();	}
		if (groupfile != "")		{		readGroup();	}
		if (listfile != "")			{		readList();		}
		if (taxfile != "")			{		readTax();		}
		if (sharedfile != "")		{		readShared();	}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
				
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string current = "";
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
			}
			
			itTypes = outputTypes.find("name");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
			}
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTaxonomyFile(current); }
			}
			
			itTypes = outputTypes.find("shared");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
			}
		}
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveGroupsCommand::readFasta(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "pick" + m->getExtension(fastafile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) == 0) {
					wroteSomething = true;
					currSeq.printSequence(out); 
				}else { 
					//if you are not in the accnos file check if you are a name that needs to be changed
					map<string, string>::iterator it = uniqueToRedundant.find(name);
					if (it != uniqueToRedundant.end()) {
						wroteSomething = true;
						currSeq.setName(it->second);
						currSeq.printSequence(out);
					}else { removedCount++; }
				}
			}
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove."); m->mothurOutEndLine();  }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your fasta file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveGroupsCommand::readShared(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sharedfile);  }
		
		//get group names from sharedfile so we can set Groups to the groupNames we want to keep
		//that way we can take advantage of the reads in inputdata and sharedRabundVector
		InputData* tempInput = new InputData(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = tempInput->getSharedRAbundVectors();
	
		//save m->Groups
		vector<string> allGroupsNames = m->getAllGroups();
		vector<string> mothurOutGroups = m->getGroups();
		
		vector<string> groupsToKeep;
		for (int i = 0; i < allGroupsNames.size(); i++) {
			if (!m->inUsersGroups(allGroupsNames[i], m->getGroups())) {
				groupsToKeep.push_back(allGroupsNames[i]);
			}
		}
		
		if (allGroupsNames.size() == groupsToKeep.size()) { m->mothurOut("Your file does not contain any groups you wish to remove."); m->mothurOutEndLine(); m->setGroups(mothurOutGroups); delete tempInput; for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  return 0; }
		
		//reset read 
		for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		delete tempInput;
		m->setGroups(groupsToKeep);
		m->clearAllGroups();
		m->names.clear();
		m->saveNextLabel = "";
		m->printedHeaders = false;
		m->currentBinLabels.clear();
		m->binLabelsInFile.clear();
		
		InputData input(sharedfile, "sharedfile");
		lookup = input.getSharedRAbundVectors();

		bool wroteSomething = false;
		
		while(lookup[0] != NULL) {
			
			string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(sharedfile)) + lookup[0]->getLabel() + ".pick" + m->getExtension(sharedfile);
			ofstream out;
			m->openOutputFile(outputFileName, out);
			outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
			
			if (m->control_pressed) { out.close();  m->mothurRemove(outputFileName);  for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; } return 0; }
			
			lookup[0]->printHeaders(out); 
			
			for (int i = 0; i < lookup.size(); i++) {
				out << lookup[i]->getLabel() << '\t' << lookup[i]->getGroup() << '\t';
				lookup[i]->print(out);
				wroteSomething = true;
				
			}			
			
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input.getSharedRAbundVectors();
			
			out.close();
		}
		
		
		m->setGroups(mothurOutGroups);
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only the groups you wish to remove."); m->mothurOutEndLine();  }
		
		string groupsString = "";
		for (int i = 0; i < Groups.size()-1; i++) {	groupsString += Groups[i] + ", "; }
		groupsString += Groups[Groups.size()-1];
		
		m->mothurOut("Removed groups: " + groupsString + " from your shared file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveGroupsCommand::readList(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + "pick" +  m->getExtension(listfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(listfile, in);
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			
			removedCount = 0;
			
			//read in list vector
			ListVector list(in);
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list.getLabel());
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
				if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
				
				//parse out names that are in accnos file
				string binnames = list.get(i);
				
				string newNames = "";
				while (binnames.find_first_of(',') != -1) { 
					string name = binnames.substr(0,binnames.find_first_of(','));
					binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
					
					//if that name is in the .accnos file, add it
					if (names.count(name) == 0) {  newNames += name + ",";  }
					else {
						//if you are not in the accnos file check if you are a name that needs to be changed
						map<string, string>::iterator it = uniqueToRedundant.find(name);
						if (it != uniqueToRedundant.end()) {
							newNames += it->second + ",";
						}else { removedCount++; }
					}
				}
				
				//get last name
				if (names.count(binnames) == 0) {  newNames += binnames + ",";  }
				else { //if you are not in the accnos file check if you are a name that needs to be changed
					map<string, string>::iterator it = uniqueToRedundant.find(binnames);
					if (it != uniqueToRedundant.end()) {
						newNames += it->second + ",";
					}else { removedCount++; }
				}
				
				//if there are names in this bin add to new list
				if (newNames != "") {  
					newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
					newList.push_back(newNames);	
				}
			}
			
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.print(out);
			}
			
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove."); m->mothurOutEndLine();  }
		outputTypes["list"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your list file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveGroupsCommand::readName(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(namefile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(namefile)) + "pick" + m->getExtension(namefile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			in >> firstCol;		m->gobble(in);		
			in >> secondCol;			
			
			vector<string> parsedNames;
			m->splitAtComma(secondCol, parsedNames);
						
			vector<string> validSecond;  validSecond.clear();
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
					validSecond.push_back(parsedNames[i]);
				}
			}
			
			removedCount += parsedNames.size()-validSecond.size();
			
			//if the name in the first column is in the set then print it and any other names in second column also in set
			if (names.count(firstCol) == 0) {
				
				wroteSomething = true;
				
				out << firstCol << '\t';
				
				//you know you have at least one valid second since first column is valid
				for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
				out << validSecond[validSecond.size()-1] << endl;
				
				//make first name in set you come to first column and then add the remaining names to second column
			}else {
				
				//you want part of this row
				if (validSecond.size() != 0) {
					
					wroteSomething = true;
					
					out << validSecond[0] << '\t';
					
					//you know you have at least one valid second since first column is valid
					for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
					out << validSecond[validSecond.size()-1] << endl;
					uniqueToRedundant[firstCol] = validSecond[0];
				}
			}
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove."); m->mothurOutEndLine();  }
		outputTypes["name"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your name file."); m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveGroupsCommand::readGroup(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "pick" + m->getExtension(groupfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << group << endl;
			}else {  removedCount++;  }
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove."); m->mothurOutEndLine();  }
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your group file."); m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveGroupsCommand::readTax(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(taxfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(taxfile)) + "pick" + m->getExtension(taxfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
			
			in >> name;				//read from first column
			in >> tax;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << tax << endl;
			}else {  //if you are not in the accnos file check if you are a name that needs to be changed
				map<string, string>::iterator it = uniqueToRedundant.find(name);
				if (it != uniqueToRedundant.end()) {
					wroteSomething = true;
					out << it->second << '\t' << tax << endl;
				}else { removedCount++; }  }
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove."); m->mothurOutEndLine();  }
		outputTypes["taxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your taxonomy file."); m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readAccnos(){
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
		
		m->setGroups(Groups);
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readAccnos");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveGroupsCommand::fillNames(){
	try {
		vector<string> seqs = groupMap->getNamesSeqs();
		
		for (int i = 0; i < seqs.size(); i++) {
			
			if (m->control_pressed) { return 0; }
			
			string group = groupMap->getGroup(seqs[i]);
			
			if (m->inUsersGroups(group, Groups)) {
				names.insert(seqs[i]);
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "fillNames");
		exit(1);
	}
}

//**********************************************************************************************************************



