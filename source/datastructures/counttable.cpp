//
//  counttable.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/26/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "counttable.h"

/************************************************************/
int CountTable::createTable(set<string>& n, map<string, string>& g, set<string>& gs) {
    try {
        int numGroups = 0;
        groups.clear();
        totalGroups.clear();
        indexGroupMap.clear();
        indexNameMap.clear();
        counts.clear();
        for (set<string>::iterator it = gs.begin(); it != gs.end(); it++) { groups.push_back(*it);  hasGroups = true; }
        numGroups = groups.size();
        totalGroups.resize(numGroups, 0);
        
		//sort groups to keep consistent with how we store the groups in groupmap
        sort(groups.begin(), groups.end());
        for (int i = 0; i < groups.size(); i++) {  indexGroupMap[groups[i]] = i; }
        m->setAllGroups(groups);
        
        uniques = 0;
        total = 0;
        bool error = false;
        //n contains treenames
        for (set<string>::iterator it = n.begin(); it != n.end(); it++) {
            
            if (m->control_pressed) { break; }
            
            string seqName = *it;
            
            vector<int> groupCounts; groupCounts.resize(numGroups, 0);
            map<string, string>::iterator itGroup = g.find(seqName);
            
            if (itGroup != g.end()) {   
                groupCounts[indexGroupMap[itGroup->second]] = 1; 
                totalGroups[indexGroupMap[itGroup->second]]++;
            }else {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
                m->mothurOut("[ERROR]: Your group file does not contain " + seqName + ". Please correct."); m->mothurOutEndLine();
            }
            
            map<string, int>::iterator it2 = indexNameMap.find(seqName);
            if (it2 == indexNameMap.end()) {
                if (hasGroups) {  counts.push_back(groupCounts);  }
                indexNameMap[seqName] = uniques;
                totals.push_back(1);
                total++;
                uniques++;
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine();
            }

        }
        if (error) { m->control_pressed = true; }
        else { //check for zero groups
            if (hasGroups) {
                for (int i = 0; i < totalGroups.size(); i++) {
                    if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
                }
            }
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "createTable");
		exit(1);
	}
}
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
int CountTable::createTable(string namefile, string groupfile, bool createGroup) {
    try {
        
        if (namefile == "") { m->mothurOut("[ERROR]: namefile cannot be blank when creating a count table.\n"); m->control_pressed = true; }
                                           
        GroupMap* groupMap;
        int numGroups = 0;
        groups.clear();
        totalGroups.clear();
        indexGroupMap.clear();
        indexNameMap.clear();
        counts.clear();
        map<int, string> originalGroupIndexes;
        
        if (groupfile != "") { 
            hasGroups = true;
            groupMap = new GroupMap(groupfile); groupMap->readMap();
            numGroups = groupMap->getNumGroups();
            groups = groupMap->getNamesOfGroups();
            totalGroups.resize(numGroups, 0);
        }else if(createGroup) {
            hasGroups = true;
            numGroups = 1;
            groups.push_back("Group1");
            totalGroups.resize(numGroups, 0);
        }
		//sort groups to keep consistent with how we store the groups in groupmap
        sort(groups.begin(), groups.end());
        for (int i = 0; i < groups.size(); i++) {  indexGroupMap[groups[i]] = i; }
        m->setAllGroups(groups);
        
        bool error = false;
        string name;
        uniques = 0;
        total = 0;
        
        
        //open input file
        ifstream in;
        m->openInputFile(namefile, in);
        
        int total = 0;
        while (!in.eof()) {
            if (m->control_pressed) { break; }
            
            string firstCol, secondCol;
            in >> firstCol; m->gobble(in); in >> secondCol; m->gobble(in);
            
            m->checkName(firstCol);
            m->checkName(secondCol);
            
            vector<string> names;
            m->splitAtChar(secondCol, names, ',');
            
            map<string, int> groupCounts;
            int thisTotal = 0;
            if (groupfile != "") {
                //set to 0
                for (int i = 0; i < groups.size(); i++) { groupCounts[groups[i]] = 0; }
                
                //get counts for each of the users groups
                for (int i = 0; i < names.size(); i++) {
                    string group = groupMap->getGroup(names[i]);
                    
                    if (group == "not found") { m->mothurOut("[ERROR]: " + names[i] + " is not in your groupfile, please correct."); m->mothurOutEndLine(); error=true; }
                    else {
                        map<string, int>::iterator it = groupCounts.find(group);
                        
                        //if not found, then this sequence is not from a group we care about
                        if (it != groupCounts.end()) {
                            it->second++;
                            thisTotal++;
                        }
                    }
                }
            }else if (createGroup) {
                groupCounts["Group1"]=0;
                for (int i = 0; i < names.size(); i++) {
                    string group = "Group1";
                    groupCounts["Group1"]++; thisTotal++;
                }
            }else { thisTotal = names.size();  }
            
            //if group info, then read it
            vector<int> thisGroupsCount; thisGroupsCount.resize(numGroups, 0);
            for (int i = 0; i < numGroups; i++) {  
                thisGroupsCount[i] = groupCounts[groups[i]]; 
                totalGroups[i] += thisGroupsCount[i]; 
            }
            
            map<string, int>::iterator it = indexNameMap.find(firstCol);
            if (it == indexNameMap.end()) {
                if (hasGroups) {  counts.push_back(thisGroupsCount);  }
                indexNameMap[firstCol] = uniques;
                totals.push_back(thisTotal);
                total += thisTotal;
                uniques++;
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + firstCol + ", sequence names must be unique. Please correct."); m->mothurOutEndLine(); 
            }
        }
        in.close();
		
        if (error) { m->control_pressed = true; }
        else { //check for zero groups
            if (hasGroups) {
                for (int i = 0; i < totalGroups.size(); i++) {
                    if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
                }
            }
        }
        if (groupfile != "") { delete groupMap; }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "createTable");
		exit(1);
	}
}
/************************************************************/
int CountTable::readTable(string file, bool readGroups, bool mothurRunning) {
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
        if ((columnHeaders.size() > 2) && readGroups) { hasGroups = true; numGroups = columnHeaders.size() - 2;  }
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
            
            if ((thisTotal == 0) && !mothurRunning) { error=true; m->mothurOut("[ERROR]: Your count table contains a sequence named " + name + " with a total=0. Please correct."); m->mothurOutEndLine();
            }
            
            //if group info, then read it
            vector<int> groupCounts; groupCounts.resize(numGroups, 0);
            if (columnHeaders.size() > 2) { //file contains groups
                if (readGroups) { //user wants to save them
                    for (int i = 0; i < numGroups; i++) {  int thisIndex = indexGroupMap[originalGroupIndexes[i]]; in >> groupCounts[thisIndex]; m->gobble(in); totalGroups[thisIndex] += groupCounts[thisIndex];  }
                }else { //read and discard
                    m->getline(in); m->gobble(in);
                }
            }
            
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
        else { //check for zero groups
            if (hasGroups) {
                for (int i = 0; i < totalGroups.size(); i++) {
                    if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
                }
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "readTable");
		exit(1);
	}
}
/************************************************************/
int CountTable::printTable(string file) {
    try {
        ofstream out;
        m->openOutputFile(file, out); 
		out << "Representative_Sequence\ttotal";
        for (int i = 0; i < groups.size(); i++) { out << '\t' << groups[i]; }
        out << endl;
        
        map<int, string> reverse; //use this to preserve order
        for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
        
        for (int i = 0; i < totals.size(); i++) {
            map<int, string>::iterator itR = reverse.find(i);
            
            if (itR != reverse.end()) { //will equal end if seqs were removed because remove just removes from indexNameMap
                out << itR->second << '\t' << totals[i];
                if (hasGroups) {
                    for (int j = 0; j < groups.size(); j++) {
                        out << '\t' << counts[i][j];
                    }
                }
                out << endl;
            }
        }
        /*for (map<string, int>::iterator itNames = indexNameMap.begin(); itNames != indexNameMap.end(); itNames++) {
            out << itNames->first << '\t' << totals[itNames->second];
            if (hasGroups) {
                
                for (int i = 0; i < groups.size(); i++) {
                    out << '\t' << counts[itNames->second][i];
                }
            }
            out << endl;
        }*/
        out.close();
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "printTable");
		exit(1);
	}
}
/************************************************************/
int CountTable::printHeaders(ofstream& out) {
    try {
		out << "Representative_Sequence\ttotal";
        for (int i = 0; i < groups.size(); i++) { out << '\t' << groups[i]; }
        out << endl;
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "printHeaders");
		exit(1);
	}
}
/************************************************************/
int CountTable::printSeq(ofstream& out, string seqName) {
    try {
		map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
        }else { 
            out << it->first << '\t' << totals[it->second];
            if (hasGroups) {
                for (int i = 0; i < groups.size(); i++) {
                    out << '\t' << counts[it->second][i];
                }
            }
            out << endl;
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "printSeq");
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
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
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
                m->mothurOut("[ERROR]: group " + groupName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
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
                m->mothurOut("[ERROR]: group " + groupName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                map<string, int>::iterator it2 = indexNameMap.find(seqName);
                if (it2 == indexNameMap.end()) {
                    //look for it in names of groups to see if the user accidently used the wrong file
                    if (m->inUsersGroups(seqName, groups)) {
                        m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                    }
                    m->mothurOut("[ERROR]: seq " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
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
//set the number of sequences for the seq for the group
int CountTable::setAbund(string seqName, string groupName, int num) {
    try {
        if (hasGroups) {
            map<string, int>::iterator it = indexGroupMap.find(groupName);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + groupName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                map<string, int>::iterator it2 = indexNameMap.find(seqName);
                if (it2 == indexNameMap.end()) {
                    //look for it in names of groups to see if the user accidently used the wrong file
                    if (m->inUsersGroups(seqName, groups)) {
                        m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                    }
                    m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
                }else { 
                    int oldCount = counts[it2->second][it->second];
                    counts[it2->second][it->second] = num;
                    totalGroups[it->second] += (num - oldCount);
                    total += (num - oldCount);
                    totals[it2->second] += (num - oldCount);
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->control_pressed = true; }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "set");
		exit(1);
	}
}
/************************************************************/
//add group
int CountTable::addGroup(string groupName) {
    try {        
        bool sanity = m->inUsersGroups(groupName, groups);
        if (sanity) { m->mothurOut("[ERROR]: " + groupName + " is already in the count table, cannot add again.\n"); m->control_pressed = true;  return 0; }
        
        groups.push_back(groupName);
        if (!hasGroups) { counts.resize(uniques);  }
        
        for (int i = 0; i < counts.size(); i++) { counts[i].push_back(0); }
        totalGroups.push_back(0);
        indexGroupMap[groupName] = groups.size()-1;
        map<string, int> originalGroupMap = indexGroupMap;
        
        //important to play well with others, :)
        sort(groups.begin(), groups.end());
        
        //fix indexGroupMap && totalGroups
        vector<int> newTotals; newTotals.resize(groups.size(), 0);
        for (int i = 0; i < groups.size(); i++) {  
            indexGroupMap[groups[i]] = i;  
            //find original spot of group[i]
            int index = originalGroupMap[groups[i]];
            newTotals[i] = totalGroups[index];
        }
        totalGroups = newTotals;
        
        //fix counts vectors
        for (int i = 0; i < counts.size(); i++) {
            vector<int> newCounts; newCounts.resize(groups.size(), 0);
            for (int j = 0; j < groups.size(); j++) {  
                //find original spot of group[i]
                int index = originalGroupMap[groups[j]];
                newCounts[j] = counts[i][index];
            }
            counts[i] = newCounts;
        }
        hasGroups = true;
        m->setAllGroups(groups);
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "addGroup");
		exit(1);
	}
}
/************************************************************/
//remove group
int CountTable::removeGroup(string groupName) {
    try {        
        if (hasGroups) {
            //save for later in case removing a group means we need to remove a seq.
            map<int, string> reverse;
            for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
            
            map<string, int>::iterator it = indexGroupMap.find(groupName);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + groupName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                int indexOfGroupToRemove = it->second;
                map<string, int> currentGroupIndex = indexGroupMap;
                vector<string> newGroups;
                for (int i = 0; i < groups.size(); i++) {
                    if (groups[i] != groupName) { 
                        newGroups.push_back(groups[i]);
                        indexGroupMap[groups[i]] = newGroups.size()-1;
                    }
                }
                indexGroupMap.erase(groupName);
                groups = newGroups;
                totalGroups.erase(totalGroups.begin()+indexOfGroupToRemove);
                
                int thisIndex = 0;
                map<string, int> newIndexNameMap;
                for (int i = 0; i < counts.size(); i++) {
                    int num = counts[i][indexOfGroupToRemove];
                    counts[i].erase(counts[i].begin()+indexOfGroupToRemove);
                    totals[i] -= num;
                    total -= num;
                    if (totals[i] == 0) { //your sequences are only from the group we want to remove, then remove you.
                        counts.erase(counts.begin()+i);
                        totals.erase(totals.begin()+i);
                        uniques--;
                        i--;
                    }
                    newIndexNameMap[reverse[thisIndex]] = i;
                    thisIndex++;
                }
                indexNameMap = newIndexNameMap;
                
                if (groups.size() == 0) { hasGroups = false; }
            }
        }else { m->mothurOut("[ERROR]: your count table does not contain group information, can not remove group " + groupName + ".\n"); m->control_pressed = true; }
    
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "removeGroup");
		exit(1);
	}
}
/************************************************************/
//vector of groups for the seq
vector<string> CountTable::getGroups(string seqName) {
    try {
        vector<string> thisGroups;
        if (hasGroups) {
            vector<int> thisCounts = getGroupCounts(seqName);
            for (int i = 0; i < thisCounts.size(); i++) {  
                if (thisCounts[i] != 0) {  thisGroups.push_back(groups[i]); }
            } 
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->control_pressed = true; }
        
        return thisGroups;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getGroups");
		exit(1);
	}
}
/************************************************************/
//total number of seqs represented by seq
int CountTable::renameSeq(string oldSeqName, string newSeqName) {
    try {
        
        map<string, int>::iterator it = indexNameMap.find(oldSeqName);
        if (it == indexNameMap.end()) {
            if (hasGroupInfo()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(oldSeqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + oldSeqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
            }
            m->mothurOut("[ERROR]: " + oldSeqName + " is not in your count table. Please correct.\n"); m->control_pressed = true;
        }else {  
            int index = it->second;
            indexNameMap.erase(it);
            indexNameMap[newSeqName] = index;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "renameSeq");
		exit(1);
	}
}

/************************************************************/
//total number of seqs represented by seq
int CountTable::getNumSeqs(string seqName) {
    try {
                
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if (hasGroupInfo()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
            }
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
//set total number of seqs represented by seq
int CountTable::setNumSeqs(string seqName, int abund) {
    try {
        
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->control_pressed = true; return -1;
        }else {
            int diff = totals[it->second] - abund;
            totals[it->second] = abund;
            total-=diff;
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
            if (hasGroupInfo()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
            }
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
        
        return 1;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "push_back");
		exit(1);
	}
}
/************************************************************/
//remove sequence
int CountTable::remove(string seqName) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it != indexNameMap.end()) {
            uniques--;
            if (hasGroups){ //remove this sequences counts from group totals
                for (int i = 0; i < totalGroups.size(); i++) {  totalGroups[i] -= counts[it->second][i];  counts[it->second][i] = 0; }
            }
            int thisTotal = totals[it->second]; totals[it->second] = 0;
            total -= thisTotal;
            indexNameMap.erase(it);
        }else {
            if (hasGroupInfo()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
            }
            m->mothurOut("[ERROR]: Your count table contains does not include " + seqName + ", cannot remove."); m->mothurOutEndLine(); m->control_pressed = true;
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
        
        return thisTotal;
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
        int thisTotal = 0;
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if ((hasGroups) && (groupCounts.size() != getNumGroups())) {  m->mothurOut("[ERROR]: Your count table has a " + toString(getNumGroups()) + " groups and " + seqName + " has " + toString(groupCounts.size()) + ", please correct."); m->mothurOutEndLine(); m->control_pressed = true;  }
            
            for (int i = 0; i < getNumGroups(); i++) {   totalGroups[i] += groupCounts[i];  thisTotal += groupCounts[i]; }
            if (hasGroups) {  counts.push_back(groupCounts);  }
            indexNameMap[seqName] = uniques;
            totals.push_back(thisTotal);
            total+= thisTotal;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct."); m->mothurOutEndLine(); m->control_pressed = true;
        }
        
        return thisTotal;
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
//returns the names of all unique sequences in file mapped to their seqCounts
map<string, int> CountTable::getNameMap() {
    try {
        map<string, int> names;
        for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) {
            names[it->first] = totals[it->second];
        }
        
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getNameMap");
		exit(1);
	}
}
/************************************************************/
//returns the names of all unique sequences in file
vector<string> CountTable::getNamesOfSeqs(string group) {
    try {
        vector<string> names;
        if (hasGroups) {
            map<string, int>::iterator it = indexGroupMap.find(group);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + group + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                for (map<string, int>::iterator it2 = indexNameMap.begin(); it2 != indexNameMap.end(); it2++) {
                    if (counts[it2->second][it->second] != 0) {  names.push_back(it2->first); }
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->control_pressed = true; }
        
        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getNamesOfSeqs");
		exit(1);
	}
}
/************************************************************/
//merges counts of seq1 and seq2, saving in seq1
int CountTable::mergeCounts(string seq1, string seq2) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seq1);
        if (it == indexNameMap.end()) {
            if (hasGroupInfo()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (m->inUsersGroups(seq1, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seq1 + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                }
            }
            m->mothurOut("[ERROR]: " + seq1 + " is not in your count table. Please correct.\n"); m->control_pressed = true;
        }else { 
            map<string, int>::iterator it2 = indexNameMap.find(seq2);
            if (it2 == indexNameMap.end()) {
                if (hasGroupInfo()) {
                    //look for it in names of groups to see if the user accidently used the wrong file
                    if (m->inUsersGroups(seq2, groups)) {
                        m->mothurOut("[WARNING]: Your group or design file contains a group named " + seq2 + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group."); m->mothurOutEndLine();
                    }
                }
                m->mothurOut("[ERROR]: " + seq2 + " is not in your count table. Please correct.\n"); m->control_pressed = true;
            }else { 
                //merge data
                for (int i = 0; i < groups.size(); i++) { counts[it->second][i] += counts[it2->second][i]; }
                totals[it->second] += totals[it2->second];
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
int CountTable::copy(CountTable* ct) {
    try {
        vector<string> thisGroups = ct->getNamesOfGroups();
        for (int i = 0; i < thisGroups.size(); i++) { addGroup(thisGroups[i]); }
        vector<string> names = ct->getNamesOfSeqs();
                                                               
        for (int i = 0; i < names.size(); i++) {
            vector<int> thisCounts = ct->getGroupCounts(names[i]);
            push_back(names[i], thisCounts);
        }
                                                               
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "copy");
		exit(1);
	}
}

/************************************************************/


