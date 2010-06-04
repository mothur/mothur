/*
 *  treemap.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treemap.h"

/************************************************************/

 TreeMap::TreeMap(string filename) {
	groupFileName = filename;
	openInputFile(filename, fileHandle);
}

/************************************************************/
 TreeMap::~TreeMap(){}

/************************************************************/
void TreeMap::readMap() {
		string seqName, seqGroup;
	
		while(fileHandle){
			fileHandle >> seqName;			//read from first column
			fileHandle >> seqGroup;			//read from second column

			namesOfSeqs.push_back(seqName);
			setNamesOfGroups(seqGroup);
					
			treemap[seqName].groupname = seqGroup;	//store data in map
			
			it2 = seqsPerGroup.find(seqGroup);
			if (it2 == seqsPerGroup.end()) { //if it's a new group
				seqsPerGroup[seqGroup] = 1;
			}else {//it's a group we already have
				seqsPerGroup[seqGroup]++;
			}

			gobble(fileHandle);
		}
		fileHandle.close();
}
/************************************************************/
void TreeMap::removeSeq(string seqName) {
	
	//erase name from namesOfSeqs
	for (int i = 0; i < namesOfSeqs.size(); i++) {
		if (namesOfSeqs[i] == seqName)  {
			namesOfSeqs.erase(namesOfSeqs.begin()+i);
			break;
		}
	}
	
	//decrement sequences in this group
	string group = treemap[seqName].groupname;
	seqsPerGroup[group]--;
	
	//remove seq from treemap
	it = treemap.find(seqName);
	treemap.erase(it);
}
/************************************************************/

int TreeMap::getNumGroups() {
			
	return namesOfGroups.size();	
		
}
/************************************************************/

int TreeMap::getNumSeqs() {
			
	return namesOfSeqs.size();	
		
}

/************************************************************/

string TreeMap::getGroup(string sequenceName) {
			
	it = treemap.find(sequenceName);
	if (it != treemap.end()) { //sequence name was in group file
		return it->second.groupname;	
	}else {
		return "not found";
	}
		
}
/************************************************************/
void TreeMap::setIndex(string seq, int index) {
	it = treemap.find(seq);
	if (it != treemap.end()) { //sequence name was in group file
		treemap[seq].vectorIndex = index;	
	}else {
		treemap[seq].vectorIndex = index;
		treemap[seq].groupname = "not found";
	}
}
/************************************************************/
int TreeMap::getIndex(string seq) {
	
	it = treemap.find(seq);
	// if it is a valid sequence name then return index
	if (it != treemap.end()) { return treemap[seq].vectorIndex; }
	// if not return error code
	else { return -1; }
	
}
/************************************************************/

void TreeMap::setNamesOfGroups(string seqGroup) {
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
			}
}
/************************************************************/
bool TreeMap::isValidGroup(string groupname) {
	try {
		for (int i = 0; i < namesOfGroups.size(); i++) {
			if (groupname == namesOfGroups[i]) { return true; }
		}
		
		return false;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "isValidGroup");
		exit(1);
	}
}
/***********************************************************************/

void TreeMap::print(ostream& output){
	try {
		
		for(it = treemap.begin(); it != treemap.end(); it++){
			output << it->first << '\t' << it->second.groupname << '\t' << it->second.vectorIndex << endl;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "print");
		exit(1);
	}
}

/************************************************************/
void TreeMap::makeSim(GroupMap* groupmap) {
	try {
		//set names of groups
		namesOfGroups = groupmap->namesOfGroups;
		
		//set names of seqs to names of groups
		namesOfSeqs = groupmap->namesOfGroups;
		
		// make map where key and value are both the group name since that what the tree.shared command wants
		for (int i = 0; i < namesOfGroups.size(); i++) {
			treemap[namesOfGroups[i]].groupname = namesOfGroups[i];
			seqsPerGroup[namesOfGroups[i]] = 1;
		}
		
		numGroups = namesOfGroups.size();
		
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "makeSim");
		exit(1);
	}
}
/************************************************************/
void TreeMap::makeSim(ListVector* list) {
	try {
		//set names of groups
		namesOfGroups.clear();
		for(int i = 0; i < list->size(); i++) {
			namesOfGroups.push_back(list->get(i));
		}
		
		//set names of seqs to names of groups
		namesOfSeqs = namesOfGroups;
		
		// make map where key and value are both the group name since that what the tree.shared command wants
		for (int i = 0; i < namesOfGroups.size(); i++) {
			treemap[namesOfGroups[i]].groupname = namesOfGroups[i];
			seqsPerGroup[namesOfGroups[i]] = 1;
		}
		
		numGroups = namesOfGroups.size();
		
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "makeSim");
		exit(1);
	}
}

/************************************************************/

