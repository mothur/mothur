//
//  counttable.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/26/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "counttable.h"


/************************************************************/
int CountTable::readTable(string file) {
    try {
        filename = file;
        ifstream in;
        m->openInputFile(filename, in);
        
        string headers = m->getline(in); m->gobble(in);
        vector<string> columnHeaders = m->splitWhiteSpace(headers);
        
        int numGroups = 0;
        groups.clear();
        totalGroups.clear();
        indexGroupMap.clear();
        indexNameMap.clear();
        counts.clear();
        map<int, string> originalGroupIndexes;
        if (columnHeaders.size() > 2) { hasGroups = true; numGroups = columnHeaders.size() - 2;  }
        for (int i = 2; i < columnHeaders.size(); i++) {  groups.push_back(columnHeaders[i]);  originalGroupIndexes[i-2] = columnHeaders[i]; totalGroups.push_back(0); }
        //sort groups to keep consistent with how we store the groups in groupmap
        sort(groups.begin(), groups.end());
        for (int i = 0; i < groups.size(); i++) {  indexGroupMap[groups[i]] = i; }
        m->setAllGroups(groups);
        
        bool error = false;
        string name;
        int thisTotal;
        uniques = 0;
        total = 0;
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            in >> name >> thisTotal; m->gobble(in);
            
            //if group info, then read it
            vector<int> groupCounts; groupCounts.resize(numGroups, 0);
            for (int i = 0; i < numGroups; i++) {  int thisIndex = indexGroupMap[originalGroupIndexes[i]]; in >> groupCounts[thisIndex]; m->gobble(in); totalGroups[thisIndex] += groupCounts[thisIndex];  }
            
            map<string, int>::iterator it = indexNameMap.find(name);
            if (it == indexNameMap.end()) {
                if (hasGroups) {  counts.push_back(groupCounts);  }
                indexNameMap[name] = uniques;
                totals.push_back(thisTotal);
                total += thisTotal;
                uniques++;
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + name + ", sequence names must be unique. Please correct."); m->mothurOutEndLine(); 
            }
        }
        in.close();
        
        if (error) { m->control_pressed = true; }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "readTable");
		exit(1);
	}
}
/************************************************************/
//group counts for a seq
vector<int> CountTable::getGroupCounts(string seqName) {
    try {
        vector<int> temp;
        if (hasGroups) {
            map<string, int>::iterator it = indexNameMap.find(seqName);
            if (it == indexNameMap.end()) {
                m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                temp = counts[it->second];
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n"); m->control_pressed = true; }
        
        return temp;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getGroupCounts");
		exit(1);
	}
}
/************************************************************/
//total number of sequences for the group
int CountTable::getGroupCount(string groupName) {
    try {
        if (hasGroups) {
            map<string, int>::iterator it = indexGroupMap.find(groupName);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + groupName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                return totalGroups[it->second];
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->control_pressed = true; }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getGroupCount");
		exit(1);
	}
}
/************************************************************/
//total number of sequences for the seq for the group
int CountTable::getGroupCount(string seqName, string groupName) {
    try {
        if (hasGroups) {
            map<string, int>::iterator it = indexGroupMap.find(groupName);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + groupName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                map<string, int>::iterator it2 = indexNameMap.find(seqName);
                if (it2 == indexNameMap.end()) {
                    m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
                }else { 
                    return counts[it2->second][it->second];
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->control_pressed = true; }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getGroupCount");
		exit(1);
	}
}
/************************************************************/
//total number of seqs represented by seq
int CountTable::getNumSeqs(string seqName) {
    try {
                
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
        }else { 
            return totals[it->second];
        }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getNumSeqs");
		exit(1);
	}
}
/************************************************************/
//returns names of seqs
vector<string> CountTable::getNamesOfSeqs() {
    try {
        vector<string> names;
        for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) {
            names.push_back(it->first);
        }
                
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getNamesOfSeqs");
		exit(1);
	}
}
/************************************************************/


