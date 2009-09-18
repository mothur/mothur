/*
 *  groupmap.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/1/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "groupmap.h"

/************************************************************/

 GroupMap::GroupMap(string filename) {
	groupFileName = filename;
	openInputFile(filename, fileHandle);
	index = 0;
}

/************************************************************/
 GroupMap::~GroupMap(){}

/************************************************************/
void GroupMap::readMap() {
		string seqName, seqGroup;
	
		while(fileHandle){
			fileHandle >> seqName;			//read from first column
			fileHandle >> seqGroup;			//read from second column
			
			setNamesOfGroups(seqGroup);
						
			groupmap[seqName] = seqGroup;	//store data in map
		
			gobble(fileHandle);
		}
		fileHandle.close();
}
/************************************************************/
int GroupMap::getNumGroups() { return namesOfGroups.size();	}
/************************************************************/

string GroupMap::getGroup(string sequenceName) {
			
	it = groupmap.find(sequenceName);
	if (it != groupmap.end()) { //sequence name was in group file
		return it->second;	
	}else {
		return "not found";
	}
}

/************************************************************/
void GroupMap::setGroup(string sequenceName, string groupN) {
	groupmap[sequenceName] = groupN;
}
/************************************************************/
void GroupMap::setNamesOfGroups(string seqGroup) {
			int i, count;
			count = 0;
			for (i=0; i<namesOfGroups.size(); i++) {
				if (namesOfGroups[i] != seqGroup) {
					count++; //you have not found this group
				}else {
					break; //you already have it
				}
			}
			if (count == namesOfGroups.size()) {
				namesOfGroups.push_back(seqGroup); //new group
				groupIndex[seqGroup] = index;
				index++;
			}
}
/************************************************************/
bool GroupMap::isValidGroup(string groupname) {
	try {
		for (int i = 0; i < namesOfGroups.size(); i++) {
			if (groupname == namesOfGroups[i]) { return true; }
		}
		
		return false;
	}
	catch(exception& e) {
		errorOut(e, "GroupMap", "isValidGroup");
		exit(1);
	}
}

/************************************************************/
vector<string> GroupMap::getNamesSeqs(){
	try {
	
		vector<string> names;
		
		for (it = groupmap.begin(); it != groupmap.end(); it++) {
			names.push_back(it->first);
		}
		
		return names;
	}
	catch(exception& e) {
		errorOut(e, "GroupMap", "getNamesSeqs");
		exit(1);
	}
}
/************************************************************/

