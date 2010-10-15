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
vector<string> SharedCommand::getValidParameters(){	
	try {
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedCommand::SharedCommand(){	
	try {
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SharedCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SharedCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

SharedCommand::SharedCommand(string o) : outputDir(o) {
	try {
		globaldata = GlobalData::getInstance();
		
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
		
		//getting output filename
		filename = globaldata->inputFileName;
		if (outputDir == "") { outputDir += m->hasPath(filename); }
		
		filename = outputDir + m->getRootName(m->getSimpleName(filename));
		filename = filename + "shared";
		outputTypes["shared"].push_back(filename);
		
		m->openOutputFile(filename, out);
		pickedGroups = false;
		
		groupMap = globaldata->gGroupmap;
		
		//if hte user has not specified any groups then use them all
		if (globaldata->Groups.size() == 0) {
			groups = groupMap->namesOfGroups;
		}else{ //they have specified groups
			groups = globaldata->Groups;
			pickedGroups = true;
		}
		
		//fill filehandles with neccessary ofstreams
		int i;
		ofstream* temp;
		for (i=0; i<groups.size(); i++) {
			temp = new ofstream;
			filehandles[groups[i]] = temp;
		}
		
		//set fileroot
		fileroot = outputDir + m->getRootName(m->getSimpleName(globaldata->getListFile()));
		
		//clears file before we start to write to it below
		for (int i=0; i<groups.size(); i++) {
			remove((fileroot + groups[i] + ".rabund").c_str());
			outputNames.push_back((fileroot + groups[i] + ".rabund"));
			outputTypes["rabund"].push_back((fileroot + groups[i] + ".rabund"));
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
		
		//lookup.clear();
		string errorOff = "no error";
		//errorOff = "";
		
		//if user provided an order file containing the order the shared file should be in read it
		if (globaldata->getOrderGroupFile() != "") { readOrderFile(); }
		
		//read in listfile
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		delete read;

		input = globaldata->ginput;
		SharedList = globaldata->gSharedList;
		string lastLabel = SharedList->getLabel();
		vector<SharedRAbundVector*> lookup; 
		
		if (m->control_pressed) { 
			delete input; delete SharedList; globaldata->ginput = NULL; globaldata->gSharedList = NULL; 
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
			out.close(); remove(filename.c_str()); 
			for (int i=0; i<groups.size(); i++) {  remove((fileroot + groups[i] + ".rabund").c_str());		}
			return 1; 
		}
				
		if ((globaldata->Groups.size() == 0) && (SharedList->getNumSeqs() != groupMap->getNumSeqs())) {  //if the user has not specified any groups and their files don't match exit with error
			m->mothurOut("Your group file contains " + toString(groupMap->getNumSeqs()) + " sequences and list file contains " + toString(SharedList->getNumSeqs()) + " sequences. Please correct."); m->mothurOutEndLine(); 
			
			out.close();
			remove(filename.c_str()); //remove blank shared file you made
			
			createMisMatchFile();
			
			//delete memory
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				delete it3->second;
			}
			delete input;
			globaldata->ginput = NULL;
			delete SharedList;
			globaldata->gSharedList = NULL;
			
			return 1; 
		}
		
		//if user has specified groups make new groupfile for them
		if (globaldata->Groups.size() != 0) { //make new group file
			string groups = "";
			for (int i = 0; i < globaldata->Groups.size(); i++) {
				groups += globaldata->Groups[i] + ".";
			}
		
			string newGroupFile = outputDir + m->getRootName(m->getSimpleName(globaldata->inputFileName)) + groups + "groups";
			ofstream outGroups;
			m->openOutputFile(newGroupFile, outGroups);
		
			vector<string> names = groupMap->getNamesSeqs();
			string groupName;
			for (int i = 0; i < names.size(); i++) {
				groupName = groupMap->getGroup(names[i]);
				if (isValidGroup(groupName, globaldata->Groups)) {
					outGroups << names[i] << '\t' << groupName << endl;
				}
			}
			outGroups.close();
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;	
	
		while((SharedList != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) { 
				delete input; delete SharedList; globaldata->ginput = NULL; globaldata->gSharedList = NULL; 
				for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
				out.close(); remove(filename.c_str()); 
				for (int i=0; i<groups.size(); i++) {  remove((fileroot + groups[i] + ".rabund").c_str());		}
				return 1; 
			}
		
			if(globaldata->allLines == 1 || globaldata->labels.count(SharedList->getLabel()) == 1){
					
					lookup = SharedList->getSharedRAbundVector();
					
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					if (pickedGroups) { //check for otus with no seqs in them
						eliminateZeroOTUS(lookup);
					}
					
					if (m->control_pressed) { 
						delete input; delete SharedList; globaldata->ginput = NULL; globaldata->gSharedList = NULL; 
						for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
						for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
						out.close(); remove(filename.c_str()); 
						for (int i=0; i<groups.size(); i++) {  remove((fileroot + groups[i] + ".rabund").c_str());		}
						return 1; 
					}
					
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
						delete input; delete SharedList; globaldata->ginput = NULL; globaldata->gSharedList = NULL; 
						for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
						for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;  }
						out.close(); remove(filename.c_str()); 
						for (int i=0; i<groups.size(); i++) {  remove((fileroot + groups[i] + ".rabund").c_str());		}
						return 1; 
					}
					
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
					delete input;  globaldata->ginput = NULL; 
					for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {  delete it3->second;   }
					out.close(); remove(filename.c_str()); 
					for (int i=0; i<groups.size(); i++) {  remove((fileroot + groups[i] + ".rabund").c_str());		}
					return 1; 
			}
			
			printSharedData(lookup); //prints info to the .shared file
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			delete SharedList;
		}
		
		globaldata->gSharedList = NULL;
		
		out.close();
		
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}

		
		//change format to shared  to speed up commands
		globaldata->setFormat("sharedfile");
		globaldata->setListFile("");
		globaldata->setGroupFile("");
		globaldata->setSharedFile(filename);
		
		
		if (m->control_pressed) { 
				delete input;  globaldata->ginput = NULL; 
				remove(filename.c_str()); 
				for (int i=0; i<groups.size(); i++) {  remove((fileroot + groups[i] + ".rabund").c_str());		}
				return 1; 
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
			
			globaldata->Groups.clear();
			
			//initialize bin values
			for (int i = 0; i < thislookup.size(); i++) {
				out << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() << '\t';
				thislookup[i]->print(out);
				
				globaldata->Groups.push_back(thislookup[i]->getGroup());
				
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
			
			globaldata->Groups.clear();
			
			//loop through ordered list and print the rabund
			for (int i = 0; i < order.size(); i++) {
				myIt = myMap.find(order[i]);
				
				if(myIt != myMap.end()) { //we found it
					out << (myIt->second)->getLabel() << '\t' << (myIt->second)->getGroup() << '\t';
					(myIt->second)->print(out);
					
					globaldata->Groups.push_back((myIt->second)->getGroup());
				
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
		string outputMisMatchName = outputDir + m->getRootName(m->getSimpleName(globaldata->inputFileName));
		
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
		m->openInputFile(globaldata->getOrderGroupFile(), in);
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


