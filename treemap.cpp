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
 TreeMap::~TreeMap(){};

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

int TreeMap::getNumGroups() {
			
	return seqsPerGroup.size();	
		
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
	treemap[seq].vectorIndex = index;
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
		cout << "Standard Error: " << e.what() << " has occurred in the TreeMap class Function isValidGroup. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeMap class function isValidGroup. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		cout << "Standard Error: " << e.what() << " has occurred in the TreeMap class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TreeMap class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/************************************************************/
