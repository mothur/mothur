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
void ParseListCommand::parse(int index) {
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
			globaldata = GlobalData::getInstance();
			int count = 1;
			
			//read in listfile
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			input = globaldata->ginput;
			list = globaldata->gSharedList;

			//read in group map info.
			groupMap = new GroupMap(globaldata->getGroupFile());
			groupMap->readMap();
			
			string seq, label;
			int i;
			//create new list vectors to fill with parsed data
			for (i=0; i<groupMap->getNumGroups(); i++) {
				groupOfLists[groupMap->namesOfGroups[i]] = new SharedListVector();
			}
			
			//parses and sets each groups listvector
			while(list != NULL){
				label = list->getLabel();
				
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(label) == 1){
				
					for(i=0; i<list->size(); i++) {
						parse(i); //parses data[i] list of sequence names
						for (it=listGroups.begin(); it != listGroups.end(); it++) {  //loop through map and set new list vectors
							seq = it->second;
							seq = seq.substr(1, seq.length()); //rips off extra comma
							groupOfLists[it->first]->push_back(seq); //sets new listvector for each group
						}
						listGroups.clear();
					}
					//prints each new list file
					for (i=0; i<groupMap->getNumGroups(); i++) {
						groupOfLists[groupMap->namesOfGroups[i]]->setLabel(label);
						groupOfLists[groupMap->namesOfGroups[i]]->print(*(filehandles[groupMap->namesOfGroups[i]]));
						groupOfLists[groupMap->namesOfGroups[i]]->clear();
					}
					
					cout << label << '\t' << count << endl;
				}
				
				list = input->getSharedListVector();
				count++;
			}
			
			//set groupmap for .shared commands
			if (globaldata->gGroupmap != NULL) { delete globaldata->gGroupmap; }
			globaldata->gGroupmap = groupMap; 
			
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
	delete list;
	delete input;
	delete read;	
}
//**********************************************************************************************************************
