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
		//groupMap = new GroupMap(globaldata->getGroupFile());
		//groupMap->readMap();
		groupMap = globaldata->gGroupmap;

		//fill filehandles with neccessary ofstreams
		int i;
		ofstream* temp;
		SharedListVector* templist;
		for (i=0; i<groupMap->getNumGroups(); i++) {
			temp = new ofstream;
			templist = new SharedListVector();
			filehandles[groupMap->namesOfGroups[i]] = temp;
			mapOfLists[groupMap->namesOfGroups[i]] = templist;
		}
		
		//set fileroot
		fileroot = getRootName(globaldata->getListFile());
		
		//clears file before we start to write to it below
		for (int i=0; i<groupMap->getNumGroups(); i++) {
			remove((fileroot + groupMap->namesOfGroups[i] + ".list").c_str());
		}
	
	}
	catch(exception& e) {
		errorOut(e, "ParseListCommand", "ParseListCommand");
		exit(1);
	}
	
}
/***********************************************************************/
void ParseListCommand::parse(int index, SharedListVector* list) {
	try {
		string member, bin, groupName;
		bin = list->get(index);
		
		while (bin.find_first_of(',') != -1) {//while you still have sequences
			member = bin.substr(0,bin.find_first_of(','));
			if ((bin.find_first_of(',')+1) <= bin.length()) {  //checks to make sure you don't have comma at end of string
				bin = bin.substr(bin.find_first_of(',')+1, bin.length());
			}
			
			groupName = groupMap->getGroup(member);
			if (groupName != "not found") {
				listGroups[groupName] = listGroups[groupName] + "," + member; //adds prefix to the correct group.
			}else {
				mothurOut("Error: Sequence '" + toString(member) + "' was not found in the group file, please correct\n");
			}
		}
		
		//save last name after comma
		groupName = groupMap->getGroup(bin);
		if (groupName != "not found") {
			listGroups[groupName] = listGroups[groupName] + "," + bin; //adds prefix to the correct group.
		}else {
			mothurOut("Error: Sequence '" + toString(bin) + "' was not found in the group file, please correct\n");
		}
	}
	catch(exception& e) {
		errorOut(e, "ParseListCommand", "parse");
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
			string lastLabel = list->getLabel();
		
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = globaldata->labels;
			set<int> userLines = globaldata->lines;
		
			//parses and sets each groups listvector
			//as long as you are not at the end of the file or done wih the lines you want
			while((list != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
								
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lastLabel) == 1){
					mothurOut(list->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					process(list);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					userLines.erase(count);
				}
				
				if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					delete list;
					list = input->getSharedListVector(lastLabel);
					
					mothurOut(list->getLabel() + "\t" + toString(count)); mothurOutEndLine();
					process(list);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
				}

				
				lastLabel = list->getLabel();			
				
				delete list;
				list = input->getSharedListVector();
				count++;
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
					needToRun = true;
				}else {
					mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
				}
			}
		
			//run last line if you need to
			if (needToRun == true)  {
				delete list;
				list = input->getSharedListVector(lastLabel);
					
				mothurOut(list->getLabel() + "\t" + toString(count)); mothurOutEndLine();
				process(list);
				delete list;
			}
			
			globaldata->gSharedList = NULL;
			//delete list vectors to fill with parsed data
			for (it2 = mapOfLists.begin(); it2 != mapOfLists.end(); it2++) {
				delete it2->second;
			}
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				delete it2->second;
			}
			
			delete input;  globaldata->ginput = NULL;
			delete read;

			
			return 0;
	}
	catch(exception& e) {
		errorOut(e, "ParseListCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

ParseListCommand::~ParseListCommand(){

			
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
					mapOfLists[it->first]->push_back(seq); //sets new listvector for each group
				}
				listGroups.clear();
			}
			//prints each new list file
			for (int i=0; i<groupMap->getNumGroups(); i++) {
				openOutputFileAppend(fileroot + groupMap->namesOfGroups[i] + ".list", *(filehandles[groupMap->namesOfGroups[i]]));
				mapOfLists[groupMap->namesOfGroups[i]]->setLabel(thisList->getLabel());
				mapOfLists[groupMap->namesOfGroups[i]]->print(*(filehandles[groupMap->namesOfGroups[i]]));
				mapOfLists[groupMap->namesOfGroups[i]]->clear();
				(*(filehandles[groupMap->namesOfGroups[i]])).close();
			}

	}
	catch(exception& e) {
		errorOut(e, "ParseListCommand", "process");
		exit(1);
	}
}
