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
	m = MothurOut::getInstance();
	groupFileName = filename;
	m->openInputFile(filename, fileHandle);
	index = 0;
}

/************************************************************/
 GroupMap::~GroupMap(){}

/************************************************************/
int GroupMap::readMap() {
		string seqName, seqGroup;
		int error = 0;

		while(fileHandle){
			fileHandle >> seqName;	m->gobble(fileHandle);		//read from first column
			fileHandle >> seqGroup;			//read from second column
			
			if (m->control_pressed) {  fileHandle.close();  return 1; }
	
			setNamesOfGroups(seqGroup);
			
			it = groupmap.find(seqName);
			
			if (it != groupmap.end()) { error = 1; m->mothurOut("Your groupfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
			else {
				groupmap[seqName] = seqGroup;	//store data in map
				seqsPerGroup[seqGroup]++;  //increment number of seqs in that group
			}
			m->gobble(fileHandle);
		}
		fileHandle.close();
		return error;
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
		seqsPerGroup[seqGroup] = 0;
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
		m->errorOut(e, "GroupMap", "isValidGroup");
		exit(1);
	}
}
/************************************************************/
int GroupMap::getNumSeqs(string group) {
	try {
		
		map<string, int>::iterator itNum;
		
		itNum = seqsPerGroup.find(group);
		
		if (itNum == seqsPerGroup.end()) { return 0; }
		
		return seqsPerGroup[group];
		
	}
	catch(exception& e) {
		m->errorOut(e, "GroupMap", "getNumSeqs");
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
		m->errorOut(e, "GroupMap", "getNamesSeqs");
		exit(1);
	}
}
/************************************************************/

