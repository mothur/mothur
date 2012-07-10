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
	m = MothurOut::getInstance();
	groupFileName = filename;
	m->openInputFile(filename, fileHandle);
}

/************************************************************/
 TreeMap::~TreeMap(){}
/************************************************************/
int TreeMap::readMap(string gf) {
    try {
        groupFileName = gf;
        m->openInputFile(gf, fileHandle);
        
        string seqName, seqGroup;
        int error = 0;

        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        
        while (!fileHandle.eof()) {
            if (m->control_pressed) { fileHandle.close();  return 1; }
            
            fileHandle.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, fileHandle.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    map<string, GroupIndex>::iterator itCheck = treemap.find(seqName);
                    if (itCheck != treemap.end()) { error = 1; m->mothurOut("[WARNING]: Your groupfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        namesOfSeqs.push_back(seqName);
                        treemap[seqName].groupname = seqGroup;	//store data in map
                        
                        it2 = seqsPerGroup.find(seqGroup);
                        if (it2 == seqsPerGroup.end()) { //if it's a new group
                            seqsPerGroup[seqGroup] = 1;
                        }else {//it's a group we already have
                            seqsPerGroup[seqGroup]++;
                        }				
                    }
                    pairDone = false; 
                } 
            }
        }
        fileHandle.close();
        
        return error;
    }
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "readMap");
		exit(1);
	}
}

/************************************************************/
int TreeMap::readMap() {
    try {
        string seqName, seqGroup;
        int error = 0;
        
        string rest = "";
        char buffer[4096];
        bool pairDone = false;
        bool columnOne = true;
        
        while (!fileHandle.eof()) {
            if (m->control_pressed) { fileHandle.close();  return 1; }
            
            fileHandle.read(buffer, 4096);
            vector<string> pieces = m->splitWhiteSpace(rest, buffer, fileHandle.gcount());
            
            for (int i = 0; i < pieces.size(); i++) {
                if (columnOne) {  seqName = pieces[i]; columnOne=false; }
                else  { seqGroup = pieces[i]; pairDone = true; columnOne=true; }
                
                if (pairDone) { 
                    setNamesOfGroups(seqGroup);
                    
                    map<string, GroupIndex>::iterator itCheck = treemap.find(seqName);
                    if (itCheck != treemap.end()) { error = 1; m->mothurOut("[WARNING]: Your groupfile contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();  }
                    else {
                        namesOfSeqs.push_back(seqName);
                        treemap[seqName].groupname = seqGroup;	//store data in map
                        
                        it2 = seqsPerGroup.find(seqGroup);
                        if (it2 == seqsPerGroup.end()) { //if it's a new group
                            seqsPerGroup[seqGroup] = 1;
                        }else {//it's a group we already have
                            seqsPerGroup[seqGroup]++;
                        }				
                    }
                    pairDone = false; 
                } 
            }
        }
        fileHandle.close();
        
        return error;
    }
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "readMap");
		exit(1);
	}
}
/************************************************************/
void TreeMap::addSeq(string seqName, string seqGroup) {
	
		namesOfSeqs.push_back(seqName);
		setNamesOfGroups(seqGroup);
					
		treemap[seqName].groupname = seqGroup;	//store data in map
			
		it2 = seqsPerGroup.find(seqGroup);
		if (it2 == seqsPerGroup.end()) { //if it's a new group
			seqsPerGroup[seqGroup] = 1;
		}else {//it's a group we already have
			seqsPerGroup[seqGroup]++;
		}
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
void TreeMap::makeSim(vector<string> ThisnamesOfGroups) {
	try {
		//set names of groups
		namesOfGroups = ThisnamesOfGroups;
		
		//set names of seqs to names of groups
		namesOfSeqs = ThisnamesOfGroups;
		
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
int TreeMap::getCopy(TreeMap& copy){
	try {
         
        namesOfGroups = copy.getNamesOfGroups();
		numGroups = copy.getNumGroups();
        namesOfSeqs = copy.namesOfSeqs;
        seqsPerGroup = copy.seqsPerGroup;
        treemap = copy.treemap;
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "getCopy");
		exit(1);
	}
}
/************************************************************/
vector<string> TreeMap::getNamesSeqs(){
	try {
        
		vector<string> names;
		
        for(it = treemap.begin(); it != treemap.end(); it++){
            names.push_back(it->first);
		}
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "getNamesSeqs");
		exit(1);
	}
}
/************************************************************/
vector<string> TreeMap::getNamesSeqs(vector<string> picked){
	try {
		
		vector<string> names;
		
		for(it = treemap.begin(); it != treemap.end(); it++){
			//if you are belong to one the the groups in the picked vector add you
			if (m->inUsersGroups(it->second.groupname, picked)) {
				names.push_back(it->first);
			}
		}
		
		return names;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeMap", "getNamesSeqs");
		exit(1);
	}
}

/************************************************************/

