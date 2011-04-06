/*
 *  removeotuscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "removeotuscommand.h"
#include "inputdata.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> RemoveOtusCommand::setParameters(){	
	try {
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(plist);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(paccnos);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveOtusCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.otus command removes otus containing sequences from a specfic group or set of groups.\n";
		helpString += "It outputs a new list file containing the otus containing sequences NOT from in the those specified groups.\n";
		helpString += "The remove.otus command parameters are accnos, group, list, label and groups. The group and list parameters are required, unless you have valid current files.\n";
		helpString += "You must also provide an accnos containing the list of groups to get or set the groups parameter to the groups you wish to select.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like.  You can separate group names with dashes.\n";
		helpString += "The label parameter allows you to specify which distance you want to process.\n";
		helpString += "The remove.otus command should be in the following format: remove.otus(accnos=yourAccnos, list=yourListFile, group=yourGroupFile, label=yourLabel).\n";
		helpString += "Example remove.otus(accnos=amazon.accnos, list=amazon.fn.list, group=amazon.groups, label=0.03).\n";
		helpString += "or remove.otus(groups=pasture, list=amazon.fn.list, amazon.groups, label=0.03).\n";
		helpString += "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "getHelpString");
		exit(1);
	}
}

///**********************************************************************************************************************
RemoveOtusCommand::RemoveOtusCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "RemoveOtusCommand");
		exit(1);
	}
}

//**********************************************************************************************************************
RemoveOtusCommand::RemoveOtusCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
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
			outputTypes["group"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			
			
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
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
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
			accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  accnosfile = ""; }	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  				
				groupfile = m->getGroupFile(); 
				if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current group file and the group parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") {  				
				listfile = m->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
			}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; m->mothurOut("You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); label=""; }	
			
			if ((accnosfile == "") && (Groups.size() == 0)) { m->mothurOut("You must provide an accnos file or specify groups using the groups parameter."); m->mothurOutEndLine(); abort = true; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "RemoveOtusCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveOtusCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		groupMap = new GroupMap(groupfile);
		groupMap->readMap();
		
		//get groups you want to remove
		if (accnosfile != "") { readAccnos(); }
		
		//make sure groups are valid
		//takes care of user setting groupNames that are invalid or setting groups=all
		SharedUtil* util = new SharedUtil();
		util->setGroups(Groups, groupMap->namesOfGroups);
		delete util;
		
		if (m->control_pressed) { delete groupMap; return 0; }
		
		//read through the list file keeping any otus that contain any sequence from the groups selected
		readListGroup();
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
		
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string current = "";
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
			}
		}
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtusCommand::readListGroup(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + "pick." + label +  m->getExtension(listfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		string GroupOutputDir = outputDir;
		if (outputDir == "") {  GroupOutputDir += m->hasPath(groupfile);  }
		string outputGroupFileName = GroupOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "pick." + label  + m->getExtension(groupfile);
		
		ofstream outGroup;
		m->openOutputFile(outputGroupFileName, outGroup);
		
		InputData* input = new InputData(listfile, "list");
		ListVector* list = input->getListVector();
		string lastLabel = list->getLabel();
		
		//using first label seen if none is provided
		if (label == "") { label = lastLabel; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		bool wroteSomething = false;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			
			if (m->control_pressed) {  delete list; delete input; out.close();  outGroup.close(); remove(outputFileName.c_str());  remove(outputGroupFileName.c_str());return 0;  }
			
			if(labels.count(list->getLabel()) == 1){
				processList(list, groupMap, out, outGroup, wroteSomething);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list; 
				
				list = input->getListVector(lastLabel);
				
				processList(list, groupMap, out, outGroup, wroteSomething);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();
			
			delete list; list = NULL;
			
			//get next line to process
			list = input->getListVector();				
		}
		
		
		if (m->control_pressed) {  if (list != NULL) { delete list; } delete input; out.close(); outGroup.close(); remove(outputFileName.c_str());  remove(outputGroupFileName.c_str()); return 0;  }
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) { delete list; }
			
			list = input->getListVector(lastLabel);
			
			processList(list, groupMap, out, outGroup, wroteSomething);
			
			delete list; list = NULL;
		}
		
		out.close();
		outGroup.close();
		
		if (wroteSomething == false) {  m->mothurOut("At distance " + label + " your file ONLY contains otus containing sequences from the groups you wish to remove."); m->mothurOutEndLine();  }
		outputTypes["list"].push_back(outputFileName); outputNames.push_back(outputFileName);
		outputTypes["group"].push_back(outputGroupFileName); outputNames.push_back(outputGroupFileName);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveOtusCommand::processList(ListVector*& list, GroupMap*& groupMap, ofstream& out, ofstream& outGroup, bool& wroteSomething){
	try {
		
		//make a new list vector
		ListVector newList;
		newList.setLabel(list->getLabel());
		
		int numOtus = 0;
		//for each bin
		for (int i = 0; i < list->getNumBins(); i++) {
			if (m->control_pressed) { return 0; }
			
			//parse out names that are in accnos file
			string binnames = list->get(i);
			
			bool removeBin = false;
			string groupFileOutput = "";
			
			//parse names
			string individual = "";
			int length = binnames.length();
			for(int j=0;j<length;j++){
				if(binnames[j] == ','){
					string group = groupMap->getGroup(individual);
					if (group == "not found") { m->mothurOut("[ERROR]: " + individual + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
					
					if (m->inUsersGroups(group, Groups)) {  removeBin = true; break; }
					groupFileOutput += individual + "\t" + group + "\n";
					individual = "";	
					
				}
				else{  individual += binnames[j];  }
			}
			
			if (!removeBin) { 
				//get last name
				string group = groupMap->getGroup(individual);
				if (group == "not found") { m->mothurOut("[ERROR]: " + individual + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
				
				if (m->inUsersGroups(group, Groups)) {  removeBin = true; }
				groupFileOutput += individual + "\t" + group + "\n";				
				
				if (!removeBin) {
					//if there are no sequences from the groups we want to remove in this bin add to new list, output to groupfile
					newList.push_back(binnames);	
					outGroup << groupFileOutput;
				}else {
					numOtus++;
				}
			}else {
				numOtus++;
			}
			
		}
		
		//print new listvector
		if (newList.getNumBins() != 0) {
			wroteSomething = true;
			newList.print(out);
		}
		
		m->mothurOut(newList.getLabel() + " - removed " + toString(numOtus) + " of the " + toString(list->getNumBins()) + " OTUs."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveOtusCommand::readAccnos(){
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
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveOtusCommand", "readAccnos");
		exit(1);
	}
}
//**********************************************************************************************************************



