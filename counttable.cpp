//
//  counttable.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/26/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "counttable.h"

/************************************************************/
bool CountTable::testGroups(string file) {
    try {
        m = MothurOut::getInstance(); hasGroups = false; total = 0;
        ifstream in;
        m->openInputFile(file, in);
    
        string headers = m->getline(in); m->gobble(in);
        vector<string> columnHeaders = m->splitWhiteSpace(headers);
        if (columnHeaders.size() > 2) { hasGroups = true;   }
        return hasGroups;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "readTable");
		exit(1);
	}
}
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
            
            in >> name; m->gobble(in); in >> thisTotal; m->gobble(in);
            if (m->debug) { m->mothurOut("[DEBUG]: " + name + '\t' + toString(thisTotal) + "\n"); }
            
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
//returns unique index for sequence like get in NameAssignment
int CountTable::get(string seqName) {
    try {
        
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
        }else { return it->second; }
        
        return -1;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "get");
		exit(1);
	}
}
/************************************************************/
//add seqeunce without group info
int CountTable::push_back(string seqName) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if (hasGroups) {  m->mothurOut("[ERROR]: Your count table has groups and I have no group information for " + seqName + "."); m->mothurOutEndLine(); m->control_pressed = true;  }
            indexNameMap[seqName] = uniques;
            totals.push_back(1);
            total++;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "push_back");
		exit(1);
	}
}
/************************************************************/
//add seqeunce without group info
int CountTable::push_back(string seqName, int thisTotal) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if (hasGroups) {  m->mothurOut("[ERROR]: Your count table has groups and I have no group information for " + seqName + "."); m->mothurOutEndLine(); m->control_pressed = true;  }
            indexNameMap[seqName] = uniques;
            totals.push_back(thisTotal);
            total+=thisTotal;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "push_back");
		exit(1);
	}
}
/************************************************************/
//add sequence with group info
int CountTable::push_back(string seqName, vector<int> groupCounts) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if ((hasGroups) && (groupCounts.size() != getNumGroups())) {  m->mothurOut("[ERROR]: Your count table has a " + toString(getNumGroups()) + " groups and " + seqName + " has " + toString(groupCounts.size()) + ", please correct."); m->mothurOutEndLine(); m->control_pressed = true;  }
            int thisTotal = 0;
            for (int i = 0; i < getNumGroups(); i++) {   totalGroups[i] += groupCounts[i];  thisTotal += groupCounts[i]; }
            indexNameMap[seqName] = uniques;
            totals.push_back(thisTotal);
            total+= thisTotal;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "push_back");
		exit(1);
	}
}

/************************************************************/
//create ListVector from uniques
ListVector CountTable::getListVector() {
    try {
        ListVector list(indexNameMap.size());
        for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) { 
            if (m->control_pressed) { break; }
            list.set(it->second, it->first); 
        }
        return list;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getListVector");
		exit(1);
	}
}

/************************************************************/
//returns the names of all unique sequences in file
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
//returns names of seqs
int CountTable::mergeCounts(string seq1, string seq2) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seq1);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seq1 + " is not in your count table. Please correct.\n"); m->control_pressed = true;
        }else { 
            map<string, int>::iterator it2 = indexNameMap.find(seq2);
            if (it2 == indexNameMap.end()) {
                m->mothurOut("[ERROR]: " + seq2 + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                //merge data
                for (int i = 0; i < groups.size(); i++) {
                    counts[it->second][i] += counts[it2->second][i];
                    counts[it2->second][i] = 0;
                }
                totals[it->second] += totals[it2->second];
                totals[it2->second] = 0;
                uniques--;
                indexNameMap.erase(it2); 
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getNamesOfSeqs");
		exit(1);
	}
}

/************************************************************/


