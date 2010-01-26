/*
 *  sharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedcommand.h"

//**********************************************************************************************************************

SharedCommand::SharedCommand(string o) : outputDir(o) {
	try {
		globaldata = GlobalData::getInstance();
		
		//getting output filename
		filename = globaldata->inputFileName;
		if (outputDir == "") { outputDir += hasPath(filename); }
		
		filename = outputDir + getRootName(getSimpleName(filename));
		filename = filename + "shared";
		
		openOutputFile(filename, out);
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
		fileroot = outputDir + getRootName(getSimpleName(globaldata->getListFile()));
		
		//clears file before we start to write to it below
		for (int i=0; i<groups.size(); i++) {
			remove((fileroot + groups[i] + ".rabund").c_str());
		}

	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SharedCommand::execute(){
	try {
		
		//lookup.clear();
		string errorOff = "no error";
		//errorOff = "";
cout << globaldata->inputFileName << endl;			
		//read in listfile
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		delete read;

		input = globaldata->ginput;
		SharedList = globaldata->gSharedList;
		string lastLabel = SharedList->getLabel();
		vector<SharedRAbundVector*> lookup; 
		
		if ((globaldata->Groups.size() == 0) && (SharedList->getNumSeqs() != groupMap->getNumSeqs())) {  //if the user has not specified any groups and their files don't match exit with error
			mothurOut("Your group file contains " + toString(groupMap->getNumSeqs()) + " sequences and list file contains " + toString(SharedList->getNumSeqs()) + " sequences. Please correct."); mothurOutEndLine(); 
			
			out.close();
			remove(filename.c_str()); //remove blank shared file you made
			
			createMisMatchFile();
			
			//delete memory
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				delete it3->second;
			}
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
		
			string newGroupFile = getRootName(globaldata->inputFileName) + groups + "groups";
			ofstream outGroups;
			openOutputFile(newGroupFile, outGroups);
		
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
		
			if(globaldata->allLines == 1 || globaldata->labels.count(SharedList->getLabel()) == 1){
			
					lookup = SharedList->getSharedRAbundVector();
					if (pickedGroups) { //check for otus with no seqs in them
						eliminateZeroOTUS(lookup);
					}
					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
				
					processedLabels.insert(SharedList->getLabel());
					userLabels.erase(SharedList->getLabel());
			}
			
			if ((anyLabelsToProcess(SharedList->getLabel(), userLabels, errorOff) == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = SharedList->getLabel();
					
					delete SharedList;
					SharedList = input->getSharedListVector(lastLabel); //get new list vector to process
					
					lookup = SharedList->getSharedRAbundVector();
					if (pickedGroups) { //check for otus with no seqs in them
						eliminateZeroOTUS(lookup);
					}
					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					
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
			if (pickedGroups) { //check for otus with no seqs in them
				eliminateZeroOTUS(lookup);
			}
			mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
			
			printSharedData(lookup); //prints info to the .shared file
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			delete SharedList;
		}
		
		globaldata->gSharedList = NULL;
		
		out.close();
		
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}

		globaldata->setSharedFile(filename);
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void SharedCommand::printSharedData(vector<SharedRAbundVector*> thislookup) {
	try {
		
		//initialize bin values
		for (int i = 0; i < thislookup.size(); i++) {
//cout << "in printData " << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() <<  endl;
			out << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() << '\t';
			thislookup[i]->print(out);
			
			RAbundVector rav = thislookup[i]->getRAbundVector();
			openOutputFileAppend(fileroot + thislookup[i]->getGroup() + ".rabund", *(filehandles[thislookup[i]->getGroup()]));
			rav.print(*(filehandles[thislookup[i]->getGroup()]));
			(*(filehandles[thislookup[i]->getGroup()])).close();
		}
 
	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "printSharedData");
		exit(1);
	}
}
//**********************************************************************************************************************
void SharedCommand::eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup) {
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
	
 
	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "eliminateZeroOTUS");
		exit(1);
	}
}
//**********************************************************************************************************************
void SharedCommand::createMisMatchFile() {
	try {
		ofstream outMisMatch;
		string outputMisMatchName = getRootName(globaldata->inputFileName);
		
		//you have sequences in your list file that are not in your group file
		if (SharedList->getNumSeqs() > groupMap->getNumSeqs()) { 
			outputMisMatchName += "missing.group";
			mothurOut("For a list of names that are in your list file and not in your group file, please refer to " + outputMisMatchName + "."); mothurOutEndLine();
			
			openOutputFile(outputMisMatchName, outMisMatch);
			
			//go through list and if group returns "not found" output it
			for (int i = 0; i < SharedList->getNumBins(); i++) {
			
				string names = SharedList->get(i); 
				
				while (names.find_first_of(',') != -1) { 
					string name = names.substr(0,names.find_first_of(','));
					names = names.substr(names.find_first_of(',')+1, names.length());
					string group = groupMap->getGroup(name);
				
					if(group == "not found") {	outMisMatch << name << endl;  }
				}
			
				//get last name
				string group = groupMap->getGroup(names);
				if(group == "not found") {	outMisMatch << names << endl;  }				
			}
			
			outMisMatch.close();
			
		
		}else {//you have sequences in your group file that are not in you list file
			
			outputMisMatchName += "missing.name";
			mothurOut("For a list of names that are in your group file and not in your list file, please refer to " + outputMisMatchName + "."); mothurOutEndLine();
			
			map<string, string> namesInList;
			
			//go through listfile and get names
			for (int i = 0; i < SharedList->getNumBins(); i++) {
				
				string names = SharedList->get(i); 
		
				while (names.find_first_of(',') != -1) { 
					string name = names.substr(0,names.find_first_of(','));
					names = names.substr(names.find_first_of(',')+1, names.length());
					
					namesInList[name] = name;
				}
				
				//get last name
				namesInList[names] = names;				
			}
			
			//get names of sequences in groupfile
			vector<string> seqNames = groupMap->getNamesSeqs();
		
			map<string, string>::iterator itMatch;
			
			openOutputFile(outputMisMatchName, outMisMatch);
			
			//loop through names in seqNames and if they aren't in namesIn list output them
			for (int i = 0; i < seqNames.size(); i++) {
				
				itMatch = namesInList.find(seqNames[i]);
				
				if (itMatch == namesInList.end()) {
				
					outMisMatch << seqNames[i] << endl; 
				}
			}		
			outMisMatch.close();
		}
 
	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "createMisMatchFile");
		exit(1);
	}
}

//**********************************************************************************************************************

SharedCommand::~SharedCommand(){
	//delete list;
	
	
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
		errorOut(e, "SharedCommand", "isValidGroup");
		exit(1);
	}
}
/************************************************************/


