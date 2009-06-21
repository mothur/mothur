/*
 *  parselistcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "parselistcommand.h"

//**********************************************************************************************************************
ParseListCommand::ParseListCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		//read in group map info.
		groupMap = new GroupMap(globaldata->getGroupFile());
		groupMap->readMap();

		//fill filehandles with neccessary ofstreams
		int i;
		ofstream* temp;
		for (i=0; i<groupMap->getNumGroups(); i++) {
			temp = new ofstream;
			filehandles[groupMap->namesOfGroups[i]] = temp;
		}
		
		//set fileroot
		fileroot = getRootName(globaldata->getListFile());
		
		//open output list files
		for (i=0; i<groupMap->getNumGroups(); i++) {//opens an output file for each group
			openOutputFile(fileroot + groupMap->namesOfGroups[i] + ".list", *(filehandles[groupMap->namesOfGroups[i]]));
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParseListCommand class Function ParseListCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParseListCommand class function ParseListCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/
void ParseListCommand::parse(int index, SharedListVector* list) {
	try {
		string prefix, suffix, groupsName;
		suffix = list->get(index);
	
		while (suffix.find_first_of(',') != -1) {//while you still have sequences
			prefix = suffix.substr(0,suffix.find_first_of(','));
			if ((suffix.find_first_of(',')+1) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
				suffix = suffix.substr(suffix.find_first_of(',')+1, suffix.length());
			}
			
			groupsName = groupMap->getGroup(prefix);
			if (groupsName != "not found") {
				listGroups[groupsName] = listGroups[groupsName] + "," + prefix; //adds prefix to the correct group.
			}else {
				cerr << "Error: Sequence '" << prefix << "' was not found in the group file, please correct\n";
			}
		}
		
		//save last name after comma
		groupsName = groupMap->getGroup(suffix);
		if (groupsName != "not found") {
			listGroups[groupsName] = listGroups[groupsName] + "," + suffix; //adds prefix to the correct group.
		}else {
			cerr << "Error: Sequence '" << suffix << "' was not found in the group file, please correct\n";
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParseListCommand class Function parse. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParseListCommand class function parse. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

int ParseListCommand::execute(){
	try{
			int count = 1;
			
			//read in listfile
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			input = globaldata->ginput;
			list = globaldata->gSharedList;
			SharedListVector* lastList = list;
		
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = globaldata->labels;
			set<int> userLines = globaldata->lines;
			
			//create new list vectors to fill with parsed data
			for (int i=0; i<groupMap->getNumGroups(); i++) {
				groupOfLists[groupMap->namesOfGroups[i]] = new SharedListVector();
			}
			
						
			//parses and sets each groups listvector
			//as long as you are not at the end of the file or done wih the lines you want
			while((list != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
								
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(list->getLabel()) == 1){
					cout << list->getLabel() << '\t' << count << endl;
					process(list);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					userLines.erase(count);
				}
				
				if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
					cout << lastList->getLabel() << '\t' << count << endl;
					process(lastList);
					
					processedLabels.insert(lastList->getLabel());
					userLabels.erase(lastList->getLabel());
				}

				if (count != 1) { delete lastList; }
				lastList = list;			

				list = input->getSharedListVector();
				count++;
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				cout << "Your file does not include the label "<< *it; 
				if (processedLabels.count(lastList->getLabel()) != 1) {
					cout << ". I will use " << lastList->getLabel() << "." << endl;
					needToRun = true;
				}else {
					cout << ". Please refer to " << lastList->getLabel() << "." << endl;
				}
			}
		
			//run last line if you need to
			if (needToRun == true)  {
				cout << lastList->getLabel() << '\t' << count << endl;
				process(lastList);
			}
			delete lastList;

			//set groupmap for .shared commands
			if (globaldata->gGroupmap != NULL) { delete globaldata->gGroupmap; }
			globaldata->gGroupmap = groupMap; 
			
			//close files
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { 
				ofstream* temp = it3->second;
				(*temp).close(); 
				delete it3->second;
			}
			
			//delete list vectors to fill with parsed data
			for (it2 = groupOfLists.begin(); it2 != groupOfLists.end(); it2++) {
				delete it2->second;
			}

			return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParseListCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParseListCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************

ParseListCommand::~ParseListCommand(){
	
		globaldata->gSharedList = NULL;
		delete input;  globaldata->ginput = NULL;
		delete read;
	
}
//**********************************************************************************************************************
void ParseListCommand::process(SharedListVector* thisList) {
	try {
			string seq;

			for(int i=0; i<thisList->size(); i++) {
				parse(i, thisList); //parses data[i] list of sequence names
				for (it=listGroups.begin(); it != listGroups.end(); it++) {  //loop through map and set new list vectors
					seq = it->second;
					seq = seq.substr(1, seq.length()); //rips off extra comma
					groupOfLists[it->first]->push_back(seq); //sets new listvector for each group
				}
				listGroups.clear();
			}
			//prints each new list file
			for (int i=0; i<groupMap->getNumGroups(); i++) {
				groupOfLists[groupMap->namesOfGroups[i]]->setLabel(thisList->getLabel());
				groupOfLists[groupMap->namesOfGroups[i]]->print(*(filehandles[groupMap->namesOfGroups[i]]));
				groupOfLists[groupMap->namesOfGroups[i]]->clear();
			}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParseListCommand class Function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParseListCommand class function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
