/*
 *  sharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedcommand.h"
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSharedRabunds(SharedRAbundVector* left, SharedRAbundVector* right){
	return (left->getGroup() < right->getGroup());	
}
//**********************************************************************************************************************
vector<string> SharedCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(plist);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pgroup);
		CommandParameter pordergroup("ordergroup", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pordergroup);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SharedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.shared command reads a list and group file and creates a shared file, as well as a rabund file for each group.\n";
		helpString += "The make.shared command parameters are list, group, ordergroup, groups and label. list and group are required unless a current file is available.\n";
		helpString += "The groups parameter allows you to indicate which groups you want to include, group names should be separated by dashes. ex. groups=A-B-C. Default is all groups in your groupfile.\n";
		helpString += "The label parameter allows you to indicate which labels you want to include, label names should be separated by dashes. Default is all labels in your list file.\n";
		helpString += "The ordergroup parameter allows you to indicate the order of the groups in the sharedfile, by default the groups are listed alphabetically.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedCommand::SharedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedCommand::SharedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
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
			 
				 it = parameters.find("ordergroup");
				 //user has given a template file
				 if(it != parameters.end()){ 
					 path = m->hasPath(it->second);
					 //if the user has not given a path then, add inputdir. else leave path alone.
					 if (path == "") {	parameters["ordergroup"] = inputDir + it->second;		}
				 }
			 }
			 
			 
			 //if the user changes the output directory command factory will send this info to us in the output parameter 
			 outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			 
			 //check for required parameters
			 listfile = validParameter.validFile(parameters, "list", true);
			 if (listfile == "not open") { listfile = ""; abort = true; }
			 else if (listfile == "not found") { 
				 listfile = m->getListFile(); 
				 if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				 else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			 }else { m->setListFile(listfile); }	
							
			 ordergroupfile = validParameter.validFile(parameters, "ordergroup", true);
			 if (ordergroupfile == "not open") { abort = true; }	
			 else if (ordergroupfile == "not found") { ordergroupfile = ""; }
			 			 
			 groupfile = validParameter.validFile(parameters, "group", true);
			 if (groupfile == "not open") { groupfile = ""; abort = true; }	
			 else if (groupfile == "not found") { 
				 groupfile = m->getGroupFile(); 
				 if (groupfile != "") { 
					 m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine();
					 groupMap = new GroupMap(groupfile);
					 
					 int error = groupMap->readMap();
					 if (error == 1) { abort = true; }
					 m->namesOfGroups = groupMap->namesOfGroups;
				 }
				 else { 	m->mothurOut("You have no current group file and the group parameter is required."); m->mothurOutEndLine(); abort = true; }
			 }else {  
				 groupMap = new GroupMap(groupfile);
			 
				 int error = groupMap->readMap();
				 if (error == 1) { abort = true; }
				 m->namesOfGroups = groupMap->namesOfGroups;
				 m->setGroupFile(groupfile);
			 }
			 
			 string groups = validParameter.validFile(parameters, "groups", false);			
			 if (groups == "not found") { groups = ""; }
			 else { 
				 m->splitAtDash(groups, Groups);
				 m->Groups = Groups;
			 }
			 
			 //check for optional parameter and set defaults
			 // ...at some point should added some additional type checking...
			 string label = validParameter.validFile(parameters, "label", false);			
			 if (label == "not found") { label = ""; }
			 else { 
				 if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				 else { allLines = 1;  }
			 }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SharedCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//getting output filename
		filename = listfile;
		
		if (outputDir == "") { outputDir += m->hasPath(filename); }
		
		filename = outputDir + m->getRootName(m->getSimpleName(filename));
		filename = filename + "shared";
		outputTypes["shared"].push_back(filename);
		
		m->openOutputFile(filename, out);
		pickedGroups = false;
		
		//if hte user has not specified any groups then use them all
		if (Groups.size() == 0) {
			Groups = groupMap->namesOfGroups; m->Groups = Groups;
		}else { pickedGroups = true; }
		
		//fill filehandles with neccessary ofstreams
		int i;
		ofstream* temp;
		for (i=0; i<Groups.size(); i++) {
			temp = new ofstream;
			filehandles[Groups[i]] = temp;
		}
		
		//set fileroot
		fileroot = outputDir + m->getRootName(m->getSimpleName(listfile));
		
		//clears file before we start to write to it below
		for (int i=0; i<Groups.size(); i++) {
			remove((fileroot + Groups[i] + ".rabund").c_str());
			outputNames.push_back((fileroot + Groups[i] + ".rabund"));
			outputTypes["rabund"].push_back((fileroot + Groups[i] + ".rabund"));
		}
		
		//lookup.clear();
		string errorOff = "no error";
		//errorOff = "";
		
		//if user provided an order file containing the order the shared file should be in read it
		if (ordergroupfile != "") { readOrderFile(); }
		
		input = new InputData(listfile, "shared");
		SharedList = input->getSharedListVector();
		string lastLabel = SharedList->getLabel();
		vector<SharedRAbundVector*> lookup; 
		
		if (m->control_pressed) { 
			delete input; delete SharedList; delete groupMap; 
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
			out.close(); remove(filename.c_str()); 
			for (int i=0; i<Groups.size(); i++) {  remove((fileroot + Groups[i] + ".rabund").c_str());		}
			return 0; 
		}
				
		if ((m->Groups.size() == 0) && (SharedList->getNumSeqs() != groupMap->getNumSeqs())) {  //if the user has not specified any groups and their files don't match exit with error
			m->mothurOut("Your group file contains " + toString(groupMap->getNumSeqs()) + " sequences and list file contains " + toString(SharedList->getNumSeqs()) + " sequences. Please correct."); m->mothurOutEndLine(); 
			
			out.close();
			remove(filename.c_str()); //remove blank shared file you made
			
			createMisMatchFile();
			
			//delete memory
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				delete it3->second;
			}
		
			delete input; delete SharedList; delete groupMap; 
			
			return 0; 
		}
		
		//if user has specified groups make new groupfile for them
		if (pickedGroups) { //make new group file
			string groups = "";
			if (m->Groups.size() < 4) {
				for (int i = 0; i < m->Groups.size(); i++) {
					groups += m->Groups[i] + ".";
				}
			}else { groups = "merge"; }
		
			string newGroupFile = outputDir + m->getRootName(m->getSimpleName(listfile)) + groups + "groups";
			outputTypes["group"].push_back(newGroupFile); 
			outputNames.push_back(newGroupFile);
			ofstream outGroups;
			m->openOutputFile(newGroupFile, outGroups);
		
			vector<string> names = groupMap->getNamesSeqs();
			string groupName;
			for (int i = 0; i < names.size(); i++) {
				groupName = groupMap->getGroup(names[i]);
				if (isValidGroup(groupName, m->Groups)) {
					outGroups << names[i] << '\t' << groupName << endl;
				}
			}
			outGroups.close();
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;	
	
		while((SharedList != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) { 
				delete input; delete SharedList; delete groupMap;
				for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
				out.close(); remove(filename.c_str()); 
				for (int i=0; i<Groups.size(); i++) {  remove((fileroot + Groups[i] + ".rabund").c_str());		}
				return 0; 
			}
		
			if(allLines == 1 || labels.count(SharedList->getLabel()) == 1){
					
					lookup = SharedList->getSharedRAbundVector();
					
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					if (pickedGroups) { //check for otus with no seqs in them
						eliminateZeroOTUS(lookup);
					}
					
					if (m->control_pressed) { 
						delete input; delete SharedList; delete groupMap; 
						for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
						for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
						out.close(); remove(filename.c_str()); 
						for (int i=0; i<Groups.size(); i++) {  remove((fileroot + Groups[i] + ".rabund").c_str());		}
						return 0; 
					}
					
					if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
				
					processedLabels.insert(SharedList->getLabel());
					userLabels.erase(SharedList->getLabel());
			}
			
			if ((m->anyLabelsToProcess(SharedList->getLabel(), userLabels, errorOff) == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = SharedList->getLabel();
					
					delete SharedList;
					SharedList = input->getSharedListVector(lastLabel); //get new list vector to process
					
					lookup = SharedList->getSharedRAbundVector();
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					if (pickedGroups) { //check for otus with no seqs in them
						eliminateZeroOTUS(lookup);
					}
					
					
					if (m->control_pressed) { 
						delete input; delete SharedList; delete groupMap; 
						for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
						for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
						out.close(); remove(filename.c_str()); 
						for (int i=0; i<Groups.size(); i++) {  remove((fileroot + Groups[i] + ".rabund").c_str());		}
						return 0; 
					}
					
					if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
					
					processedLabels.insert(SharedList->getLabel());
					userLabels.erase(SharedList->getLabel());
					
					//restore real lastlabel to save below
					SharedList->setLabel(saveLabel);
			}
			
		
			lastLabel = SharedList->getLabel();
				
			delete SharedList;
			SharedList = input->getSharedListVector(); //get new list vector to process
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			if (processedLabels.count(lastLabel) != 1) {
				needToRun = true;
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (SharedList != NULL) {	delete SharedList;	}
			SharedList = input->getSharedListVector(lastLabel); //get new list vector to process
					
			lookup = SharedList->getSharedRAbundVector();
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			if (pickedGroups) { //check for otus with no seqs in them
				eliminateZeroOTUS(lookup);
			}
			
			if (m->control_pressed) { 
				delete input;  delete groupMap;
					for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;   }
					out.close(); remove(filename.c_str()); 
					for (int i=0; i<Groups.size(); i++) {  remove((fileroot + Groups[i] + ".rabund").c_str());		}
					return 0; 
			}
			
			if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
			printSharedData(lookup); //prints info to the .shared file
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			delete SharedList;
		}
		
		out.close();
		
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}

		delete input; delete groupMap;
		
		if (m->control_pressed) { 
				remove(filename.c_str()); 
				for (int i=0; i<Groups.size(); i++) {  remove((fileroot + Groups[i] + ".rabund").c_str());		}
				return 0; 
		}
		
		//set rabund file as new current rabundfile
		string current = "";
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}	
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOut(filename); m->mothurOutEndLine();
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void SharedCommand::printSharedData(vector<SharedRAbundVector*> thislookup) {
	try {
		
		if (order.size() == 0) { //user has not specified an order so do aplabetically
			sort(thislookup.begin(), thislookup.end(), compareSharedRabunds);
			
			m->Groups.clear();
			
			//initialize bin values
			for (int i = 0; i < thislookup.size(); i++) {
				out << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() << '\t';
				thislookup[i]->print(out);
				
				m->Groups.push_back(thislookup[i]->getGroup());
				
				RAbundVector rav = thislookup[i]->getRAbundVector();
				m->openOutputFileAppend(fileroot + thislookup[i]->getGroup() + ".rabund", *(filehandles[thislookup[i]->getGroup()]));
				rav.print(*(filehandles[thislookup[i]->getGroup()]));
				(*(filehandles[thislookup[i]->getGroup()])).close();
			}
		}else{
			//create a map from groupName to each sharedrabund
			map<string, SharedRAbundVector*> myMap;
			map<string, SharedRAbundVector*>::iterator myIt;
			
			for (int i = 0; i < thislookup.size(); i++) {
				myMap[thislookup[i]->getGroup()] = thislookup[i];
			}
			
			m->Groups.clear();
			
			//loop through ordered list and print the rabund
			for (int i = 0; i < order.size(); i++) {
				myIt = myMap.find(order[i]);
				
				if(myIt != myMap.end()) { //we found it
					out << (myIt->second)->getLabel() << '\t' << (myIt->second)->getGroup() << '\t';
					(myIt->second)->print(out);
					
					m->Groups.push_back((myIt->second)->getGroup());
				
					RAbundVector rav = (myIt->second)->getRAbundVector();
					m->openOutputFileAppend(fileroot + (myIt->second)->getGroup() + ".rabund", *(filehandles[(myIt->second)->getGroup()]));
					rav.print(*(filehandles[(myIt->second)->getGroup()]));
					(*(filehandles[(myIt->second)->getGroup()])).close();
				}else{
					m->mothurOut("Can't find shared info for " + order[i] + ", skipping."); m->mothurOutEndLine();
				}
			}
		
		}
 
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "printSharedData");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedCommand::eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup) {
	try {
		
		vector<SharedRAbundVector*> newLookup;
		for (int i = 0; i < thislookup.size(); i++) {
			SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(thislookup[i]->getLabel());
			temp->setGroup(thislookup[i]->getGroup());
			newLookup.push_back(temp);
		}
		
		//for each bin
		for (int i = 0; i < thislookup[0]->getNumBins(); i++) {
			if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
		
			//look at each sharedRabund and make sure they are not all zero
			bool allZero = true;
			for (int j = 0; j < thislookup.size(); j++) {
				if (thislookup[j]->getAbundance(i) != 0) { allZero = false;  break;  }
			}
			
			//if they are not all zero add this bin
			if (!allZero) {
				for (int j = 0; j < thislookup.size(); j++) {
					newLookup[j]->push_back(thislookup[j]->getAbundance(i), thislookup[j]->getGroup());
				}
			}
			//else{  cout << "bin # " << i << " is all zeros" << endl;  }
		}
	
		for (int j = 0; j < thislookup.size(); j++) {  delete thislookup[j];  }
		thislookup = newLookup;
		
		return 0;
 
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "eliminateZeroOTUS");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedCommand::createMisMatchFile() {
	try {
		ofstream outMisMatch;
		string outputMisMatchName = outputDir + m->getRootName(m->getSimpleName(listfile));
		
		//you have sequences in your list file that are not in your group file
		if (SharedList->getNumSeqs() > groupMap->getNumSeqs()) { 
			outputMisMatchName += "missing.group";
			m->mothurOut("For a list of names that are in your list file and not in your group file, please refer to " + outputMisMatchName + "."); m->mothurOutEndLine();
			
			m->openOutputFile(outputMisMatchName, outMisMatch);
			
			map<string, string> listNames;
			map<string, string>::iterator itList;
			
			//go through list and if group returns "not found" output it
			for (int i = 0; i < SharedList->getNumBins(); i++) {
				if (m->control_pressed) { outMisMatch.close(); remove(outputMisMatchName.c_str()); return 0; } 
			
				string names = SharedList->get(i); 
				
				while (names.find_first_of(',') != -1) { 
					string name = names.substr(0,names.find_first_of(','));
					names = names.substr(names.find_first_of(',')+1, names.length());
					string group = groupMap->getGroup(name);
					
					if(group == "not found") {	outMisMatch << name << endl;  }
					
					itList = listNames.find(name);
					if (itList != listNames.end()) {  m->mothurOut(name + " is in your list file more than once.  Sequence names must be unique. please correct."); m->mothurOutEndLine(); }
					else { listNames[name] = name; }
				}
			
				//get last name
				string group = groupMap->getGroup(names);
				if(group == "not found") {	outMisMatch << names << endl;  }	
				
				itList = listNames.find(names);
				if (itList != listNames.end()) {  m->mothurOut(names + " is in your list file more than once.  Sequence names must be unique. please correct."); m->mothurOutEndLine(); }
				else { listNames[names] = names; }

			}
			
			outMisMatch.close();
			
		
		}else {//you have sequences in your group file that are not in you list file
			
			outputMisMatchName += "missing.name";
			m->mothurOut("For a list of names that are in your group file and not in your list file, please refer to " + outputMisMatchName + "."); m->mothurOutEndLine();
			
			map<string, string> namesInList;
			map<string, string>::iterator itList;
			
			//go through listfile and get names
			for (int i = 0; i < SharedList->getNumBins(); i++) {
				if (m->control_pressed) {  return 0; } 

				
				string names = SharedList->get(i); 
		
				while (names.find_first_of(',') != -1) { 
					string name = names.substr(0,names.find_first_of(','));
					names = names.substr(names.find_first_of(',')+1, names.length());
					
					itList = namesInList.find(name);
					if (itList != namesInList.end()) {  m->mothurOut(name + " is in your list file more than once.  Sequence names must be unique. please correct."); m->mothurOutEndLine(); }

					namesInList[name] = name;
					
				}
				
				itList = namesInList.find(names);
				if (itList != namesInList.end()) {  m->mothurOut(names + " is in your list file more than once.  Sequence names must be unique. please correct."); m->mothurOutEndLine(); }

				//get last name
				namesInList[names] = names;				
			}
			
			//get names of sequences in groupfile
			vector<string> seqNames = groupMap->getNamesSeqs();
		
			map<string, string>::iterator itMatch;
			
			m->openOutputFile(outputMisMatchName, outMisMatch);
			
			//loop through names in seqNames and if they aren't in namesIn list output them
			for (int i = 0; i < seqNames.size(); i++) {
				if (m->control_pressed) { outMisMatch.close(); remove(outputMisMatchName.c_str()); return 0; } 
				
				itMatch = namesInList.find(seqNames[i]);
				
				if (itMatch == namesInList.end()) {
				
					outMisMatch << seqNames[i] << endl; 
				}
			}		
			outMisMatch.close();
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "createMisMatchFile");
		exit(1);
	}
}

//**********************************************************************************************************************

SharedCommand::~SharedCommand(){
	//delete list;
	
	
}
//**********************************************************************************************************************
int SharedCommand::readOrderFile() {
	try {
		//remove old names
		order.clear();
		
		ifstream in;
		m->openInputFile(ordergroupfile, in);
		string thisGroup;
		
		while(!in.eof()){
			in >> thisGroup; m->gobble(in);
						
			order.push_back(thisGroup);
			
			if (m->control_pressed) { order.clear(); break; }
		}
		in.close();		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "readOrderFile");
		exit(1);
	}
}
//**********************************************************************************************************************

bool SharedCommand::isValidGroup(string groupname, vector<string> groups) {
	try {
		for (int i = 0; i < groups.size(); i++) {
			if (groupname == groups[i]) { return true; }
		}
		
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "isValidGroup");
		exit(1);
	}
}
/************************************************************/


