/*
 *  removerarecommand.cpp
 *  mothur
 *
 *  Created by westcott on 1/21/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "removerarecommand.h"
#include "sequence.hpp"
#include "groupmap.h"
#include "sharedutilities.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> RemoveRareCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pshared);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pnseqs("nseqs", "Number", "", "0", "", "", "",false,true); parameters.push_back(pnseqs);
		CommandParameter pbygroup("bygroup", "Boolean", "", "f", "", "", "",false,true); parameters.push_back(pbygroup);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveRareCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.rare command parameters are list, rabund, sabund, shared, group, label, groups, bygroup and nseqs.\n";
		helpString += "The remove.rare command reads one of the following file types: list, rabund, sabund or shared file. It outputs a new file after removing the rare otus.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like analyzed.  Default=all. You may separate group names with dashes.\n";
		helpString += "The label parameter is used to analyze specific labels in your input. default=all. You may separate label names with dashes.\n";
		helpString += "The bygroup parameter is only valid with the shared file. default=f, meaning remove any OTU that has nseqs or fewer sequences across all groups.\n";
		helpString += "bygroups=T means remove any OTU that has nseqs or fewer sequences in each group (if groupA has 1 sequence and group B has 100 sequences in OTUZ and nseqs=1, then set the groupA count for OTUZ to 0 and keep groupB's count at 100.) \n";
		helpString += "The nseqs parameter allows you to set the cutoff for an otu to be deemed rare. It is required.\n";
		helpString += "The remove.rare command should be in the following format: remove.rare(shared=yourSharedFile, nseqs=yourRareCutoff).\n";
		helpString += "Example remove.rare(shared=amazon.fn.shared, nseqs=2).\n";
		helpString += "Note: No spaces between parameter labels (i.e. shared), '=' and parameters (i.e.yourSharedFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveRareCommand::RemoveRareCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "RemoveRareCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveRareCommand::RemoveRareCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
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
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["shared"] = tempOutNames;	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
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
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			
			//check for file parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			else { m->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { abort = true; }
			else if (sabundfile == "not found") {  sabundfile = "";  }	
			else { m->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { abort = true; }
			else if (rabundfile == "not found") {  rabundfile = "";  }				
			else { m->setRabundFile(rabundfile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { m->setGroupFile(groupfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = "";  abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
			else { m->setSharedFile(sharedfile); }
			
			if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					listfile = m->getListFile(); 
					if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else { 
						rabundfile = m->getRabundFile(); 
						if (rabundfile != "") {  m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
						else { 
							sabundfile = m->getSabundFile(); 
							if (sabundfile != "") {  m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
							else { 
								m->mothurOut("No valid current files. You must provide a list, sabund, rabund or shared file."); m->mothurOutEndLine(); 
								abort = true;
							}
						}
					}
				}
			} 
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "all"; }
			m->splitAtDash(groups, Groups);
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			string temp = validParameter.validFile(parameters, "nseqs", false);	 
			if (temp == "not found") { m->mothurOut("nseqs is a required parameter."); m->mothurOutEndLine(); abort = true; }
			else { convert(temp, nseqs); }
			
			temp = validParameter.validFile(parameters, "bygroup", false);	 if (temp == "not found") { temp = "f"; }
			byGroup = m->isTrue(temp);
			
			if (byGroup && (sharedfile == "")) { m->mothurOut("The byGroup parameter is only valid with a shared file."); m->mothurOutEndLine(); }
			
			if ((groupfile != "") && (listfile == "")) { m->mothurOut("A groupfile is only valid with a list file."); m->mothurOutEndLine(); groupfile = ""; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "RemoveRareCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveRareCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if (m->control_pressed) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (sabundfile != "")		{		processSabund();	}
		if (rabundfile != "")		{		processRabund();	}
		if (listfile != "")			{		processList();		}
		if (sharedfile != "")		{		processShared();	}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
			
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set rabund file as new current rabundfile
			string current = "";
			itTypes = outputTypes.find("rabund");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
			}
			
			itTypes = outputTypes.find("sabund");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSabundFile(current); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
			}
			
			itTypes = outputTypes.find("shared");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
			}
		}
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveRareCommand::processList(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + "pick" +  m->getExtension(listfile);
		string outputGroupFileName = thisOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "pick" +  m->getExtension(groupfile);
		
		ofstream out, outGroup;
		m->openOutputFile(outputFileName, out);
		
		bool wroteSomething = false;
		
		//you must provide a label because the names in the listfile need to be consistent
		string thisLabel = "";
		if (allLines) { m->mothurOut("For the listfile you must select one label, using first label in your listfile."); m->mothurOutEndLine(); }
		else if (labels.size() > 1) { m->mothurOut("For the listfile you must select one label, using " + (*labels.begin()) + "."); m->mothurOutEndLine(); thisLabel = *labels.begin(); }
		else { thisLabel = *labels.begin(); }
		
		InputData input(listfile, "list");
		ListVector* list = input.getListVector();
		
		//get first one or the one we want
		if (thisLabel != "") { 	
			//use smart distancing
			set<string> userLabels; userLabels.insert(thisLabel);
			set<string> processedLabels;
			string lastLabel = list->getLabel();
			while((list != NULL) && (userLabels.size() != 0)) {
				if(userLabels.count(list->getLabel()) == 1){
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					break;
				}
				
				if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					delete list;
					list = input.getListVector(lastLabel);
					break;
				}
				lastLabel = list->getLabel();
				delete list;
				list = input.getListVector();
			}
			if (userLabels.size() != 0) { 
				m->mothurOut("Your file does not include the label " + thisLabel + ". I will use " + lastLabel + ".");  m->mothurOutEndLine();
				list = input.getListVector(lastLabel); 
			}
		}
		
		//if groupfile is given then use it
		GroupMap* groupMap;
		if (groupfile != "") { 
			groupMap = new GroupMap(groupfile); groupMap->readMap(); 
			SharedUtil util;
			util.setGroups(Groups, groupMap->namesOfGroups);
			m->openOutputFile(outputGroupFileName, outGroup);
		}
		
		
		if (list != NULL) {	
			//make a new list vector
			ListVector newList;
			newList.setLabel(list->getLabel());
			
			//for each bin
			for (int i = 0; i < list->getNumBins(); i++) {
				if (m->control_pressed) {  if (groupfile != "") { delete groupMap; outGroup.close(); m->mothurRemove(outputGroupFileName); } out.close();  m->mothurRemove(outputFileName);  return 0; }
				
				//parse out names that are in accnos file
				string binnames = list->get(i);
				vector<string> names;
				string saveBinNames = binnames;
				m->splitAtComma(binnames, names);
				
				vector<string> newGroupFile;
				if (groupfile != "") {
					vector<string> newNames;
					saveBinNames = "";
					for(int k = 0; k < names.size(); k++) {
						string group = groupMap->getGroup(names[k]);
						
						if (m->inUsersGroups(group, Groups)) {
							newGroupFile.push_back(names[k] + "\t" + group); 
								
							newNames.push_back(names[k]);	
							saveBinNames += names[k] + ",";
						}
					}
					names = newNames;
					saveBinNames = saveBinNames.substr(0, saveBinNames.length()-1);
				}

				if (names.size() > nseqs) { //keep bin
					newList.push_back(saveBinNames);
					for(int k = 0; k < newGroupFile.size(); k++) { outGroup << newGroupFile[k] << endl; }
				}
			}
			
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.print(out);
			}
		}	
		
		out.close();
		if (groupfile != "") { outGroup.close(); outputTypes["group"].push_back(outputGroupFileName); outputNames.push_back(outputGroupFileName); }
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only rare sequences."); m->mothurOutEndLine();  }
		outputTypes["list"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processSabund(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sabundfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(sabundfile)) + "pick" +  m->getExtension(sabundfile);
		outputTypes["sabund"].push_back(outputFileName); outputNames.push_back(outputFileName);

		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		InputData input(sabundfile, "sabund");
		SAbundVector* sabund = input.getSAbundVector();
		string lastLabel = sabund->getLabel();
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { delete sabund; out.close(); return 0; }
			
			if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
				
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				if (sabund->getMaxRank() > nseqs) {
					for(int i = 1; i <=nseqs; i++) {  sabund->set(i, 0); }
				}else {	sabund->clear(); }
				
				if (sabund->getNumBins() > 0) { sabund->print(out); }
			}
			
			if ((m->anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = sabund->getLabel();
				
				delete sabund;
				sabund = input.getSAbundVector(lastLabel);
				
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				if (sabund->getMaxRank() > nseqs) {
					for(int i = 1; i <=nseqs; i++) {  sabund->set(i, 0); }
				}else {	sabund->clear(); }
				
				if (sabund->getNumBins() > 0) { sabund->print(out); }
								
				//restore real lastlabel to save below
				sabund->setLabel(saveLabel);
			}		
			
			lastLabel = sabund->getLabel();			
			
			delete sabund;
			sabund = input.getSAbundVector();
		}
		
		if (m->control_pressed) {  out.close(); return 0; }	
		
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
			if (sabund != NULL) {	delete sabund;	}
			sabund = input.getSAbundVector(lastLabel);
			
			m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
			
			if (sabund->getMaxRank() > nseqs) {
				for(int i = 1; i <=nseqs; i++) {  sabund->set(i, 0); }
			}else {	sabund->clear(); }
			
			if (sabund->getNumBins() > 0) { sabund->print(out); }
			
			delete sabund;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processSabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processRabund(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(rabundfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(rabundfile)) + "pick" +  m->getExtension(rabundfile);
		outputTypes["rabund"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		InputData input(rabundfile, "rabund");
		RAbundVector* rabund = input.getRAbundVector();
		string lastLabel = rabund->getLabel();
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { delete rabund; out.close(); return 0; }
			
			if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
				
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
				
				RAbundVector newRabund; newRabund.setLabel(rabund->getLabel());
				for (int i = 0; i < rabund->getNumBins(); i++) {
					if (rabund->get(i) > nseqs) {
						newRabund.push_back(rabund->get(i));
					}
				}
				if (newRabund.getNumBins() > 0) { newRabund.print(out); }
			}
			
			if ((m->anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = rabund->getLabel();
				
				delete rabund;
				rabund = input.getRAbundVector(lastLabel);
				
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
				
				RAbundVector newRabund; newRabund.setLabel(rabund->getLabel());
				for (int i = 0; i < rabund->getNumBins(); i++) {
					if (rabund->get(i) > nseqs) {
						newRabund.push_back(rabund->get(i));
					}
				}
				if (newRabund.getNumBins() > 0) { newRabund.print(out); }				
				
				//restore real lastlabel to save below
				rabund->setLabel(saveLabel);
			}		
			
			lastLabel = rabund->getLabel();			
			
			delete rabund;
			rabund = input.getRAbundVector();
		}
		
		if (m->control_pressed) {  out.close(); return 0; }	
		
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
			if (rabund != NULL) {	delete rabund;	}
			rabund = input.getRAbundVector(lastLabel);
			
			m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
			
			RAbundVector newRabund; newRabund.setLabel(rabund->getLabel());
			for (int i = 0; i < rabund->getNumBins(); i++) {
				if (rabund->get(i) > nseqs) {
					newRabund.push_back(rabund->get(i));
				}
			}
			if (newRabund.getNumBins() > 0) { newRabund.print(out); }	
			
			delete rabund;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processRabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processShared(){
	try {
		m->Groups = Groups;
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sharedfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(sharedfile)) + "pick" +  m->getExtension(sharedfile);
		outputTypes["shared"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		InputData input(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  out.close(); return 0; }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
				processLookup(lookup, out);
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
				processLookup(lookup, out);			
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}		
			
			lastLabel = lookup[0]->getLabel();			
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
			lookup = input.getSharedRAbundVectors();
		}
		
		if (m->control_pressed) {  out.close(); return 0; }	
		
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
			for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i];	}  }
			lookup = input.getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			
			if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
			processLookup(lookup, out);	
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processSabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processLookup(vector<SharedRAbundVector*>& lookup, ofstream& out){
	try {
		
		vector<SharedRAbundVector> newRabunds;  newRabunds.resize(lookup.size());
		for (int i = 0; i < lookup.size(); i++) {  
			newRabunds[i].setGroup(lookup[i]->getGroup());
			newRabunds[i].setLabel(lookup[i]->getLabel());
		}
		
		if (byGroup) {
			
			//for each otu
			for (int i = 0; i < lookup[0]->getNumBins(); i++) {
				bool allZero = true;
				
				if (m->control_pressed) { return 0; }
				
				//for each group
				for (int j = 0; j < lookup.size(); j++) {
					
					//are you rare?
					if (lookup[j]->getAbundance(i) > nseqs) {
						newRabunds[j].push_back(lookup[j]->getAbundance(i), newRabunds[j].getGroup());
						allZero = false;
					}else {
						newRabunds[j].push_back(0, newRabunds[j].getGroup());
					}
				}
				
				//eliminates zero otus
				if (allZero) { for (int j = 0; j < newRabunds.size(); j++) {  newRabunds[j].pop_back(); } }
			}
		}else {
			//for each otu
			for (int i = 0; i < lookup[0]->getNumBins(); i++) {
				
				if (m->control_pressed) { return 0; }
				
				int totalAbund = 0;
				//get total otu abundance
				for (int j = 0; j < lookup.size(); j++) {
					newRabunds[j].push_back(lookup[j]->getAbundance(i), newRabunds[j].getGroup());
					totalAbund += lookup[j]->getAbundance(i);
				}
				
				//eliminates otus below rare cutoff
				if (totalAbund <= nseqs) { for (int j = 0; j < newRabunds.size(); j++) {  newRabunds[j].pop_back(); } }
			}
		}
		
		//do we have any otus above the rare cutoff
		if (newRabunds[0].getNumBins() != 0) { 
			for (int j = 0; j < newRabunds.size(); j++) { 
				out << newRabunds[j].getLabel() << '\t' << newRabunds[j].getGroup() << '\t';
				newRabunds[j].print(out); 
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processLookup");
		exit(1);
	}
}
//**********************************************************************************************************************




