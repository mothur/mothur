//
//  counttable.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/26/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "counttable.h"

/************************************************************/
//used by tree commands
int CountTable::createTable(set<string>& n, map<string, string>& g, set<string>& gs) {
    try {
        hasGroups = false;
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

        uniques = 0;
        total = 0;
        bool error = false;
        //n contains treenames
        for (set<string>::iterator it = n.begin(); it != n.end(); it++) {

            if (m->getControl_pressed()) { break; }

            string seqName = *it;

            vector<countTableItem> groupCounts;
            map<string, string>::iterator itGroup = g.find(seqName);

            if (itGroup != g.end()) {
                groupCounts.push_back(countTableItem(1, indexGroupMap[itGroup->second]));
                totalGroups[indexGroupMap[itGroup->second]]++;
            }else {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (util.inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n");
                }
                m->mothurOut("[ERROR]: Your group file does not contain " + seqName + ". Please correct.\n");
            }

            map<string, int>::iterator it2 = indexNameMap.find(seqName);
            if (it2 == indexNameMap.end()) {
                if (hasGroups) { counts.push_back(groupCounts); }
                indexNameMap[seqName] = uniques;
                totals.push_back(1);
                total++;
                uniques++;
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct.\n");
            }

        }
        if (error) { m->setControl_pressed(true); }
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
        vector<string> nothing;
        return testGroups(file, nothing);
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "testGroups");
        exit(1);
    }
}

/************************************************************/
bool CountTable::testGroups(string file, vector<string>& groups) {
    try {
        m = MothurOut::getInstance(); hasGroups = false; total = 0;
        ifstream in;
        util.openInputFile(file, in);

        string headers = util.getline(in); util.gobble(in);
        
        if (headers[0] == '#') { //is this a count file in compressed form
            isCompressed = true;
            
            //read headers
            headers = util.getline(in); util.gobble(in); //gets compressed group name map line
            headers = util.getline(in); util.gobble(in);
        }
        
        vector<string> columnHeaders = util.splitWhiteSpace(headers);

        if (columnHeaders.size() > 2) {
            hasGroups = true;

            for (int i = 2; i < columnHeaders.size(); i++) {
                groups.push_back(columnHeaders[i]);
            }
            //sort groups to keep consistent with how we store the groups in groupmap
            sort(groups.begin(), groups.end());
        }

        return hasGroups;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "testGroups");
		exit(1);
	}
}

/************************************************************/

bool CountTable::setNamesOfGroups(vector<string> mygroups) {
    try {
        //remove groups from table not in new groups we are setting
        for (int i = 0; i < groups.size();) {
            if (util.inUsersGroups(groups[i], mygroups)) { ++i; }
            else { removeGroup(groups[i]);  }
        }

        //add any new groups in new groups list to table
        for (int i = 0; i < mygroups.size(); i++) {
            if (util.inUsersGroups(mygroups[i], groups)) {}
            else { addGroup(mygroups[i]);  }
        }

        //false if error
        return (!m->getControl_pressed());
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "setNamesOfGroups");
        exit(1);
    }
}

/************************************************************/

int CountTable::createTable(string namefile, string groupfile, vector<string> selectedGroups, bool createGroup) {
    try {
        
        GroupMap* groupMap;
        int numGroups = 0;
        groups.clear();
        totalGroups.clear();
        indexGroupMap.clear();
        indexNameMap.clear();
        counts.clear();
        map<int, string> originalGroupIndexes;
        uniques = 0;
        total = 0;
        bool error = false;
        bool pickedGroups = false;
        if (selectedGroups.size() != 0) { pickedGroups = true; }
        
        if (groupfile != "") {
            hasGroups = true;
            groupMap = new GroupMap(groupfile); groupMap->readMap(selectedGroups);
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

        if ((namefile == "") && (groupfile == "")) { m->mothurOut("[ERROR]: No name or group file given. You must provide a name or group file to create a count file, please correct.\n");  m->setControl_pressed(true); return 0; }
        
        else if (namefile != "") {
        
            ifstream in; util.openInputFile(namefile, in);
            
            while (!in.eof()) {
                if (m->getControl_pressed()) { break; }
                
                string firstCol, secondCol;
                in >> firstCol; util.gobble(in); in >> secondCol; util.gobble(in);
                
                util.checkName(firstCol);
                util.checkName(secondCol);
                
                vector<string> names;
                util.splitAtChar(secondCol, names, ',');
                
                map<string, int> groupCounts;
                for (int i = 0; i < groups.size(); i++) { groupCounts[groups[i]] = 0;  } //initialize groupCounts
                
                int thisTotal = 0;
                if (groupfile != "") {
                    
                    //get counts for each of the users groups
                    for (int i = 0; i < names.size(); i++) {
                        string group = groupMap->getGroup(names[i]);
                        
                        if (group == "not found") {
                            if (!pickedGroups) { m->mothurOut("[ERROR]: " + names[i] + " is not in your groupfile, please correct.\n");  error=true; }
                            //else - ignore because we assume this read is from a group we are not interested in
                        }else { //this is a read from a group we want to save
                            map<string, int>::iterator it = groupCounts.find(group);
                            
                            //if not found, then this sequence is not from a group we care about
                            if (it != groupCounts.end()) { it->second++; }
                            thisTotal++;
                        }
                    }
                }else if (createGroup) {
                    thisTotal = names.size();
                    groupCounts["Group1"] = thisTotal;
                }else { thisTotal = names.size();  }
                
                //if group info, then read it
                vector<countTableItem> thisGroupsCount;
                for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) {
                    int groupIndex = indexGroupMap[it->first];
                    int abund = it->second;
                    if (abund != 0) {
                        countTableItem thisAbund(it->second, groupIndex);
                        thisGroupsCount.push_back(thisAbund);
                        totalGroups[groupIndex] += abund;
                    }
                    
                }
                
                map<string, int>::iterator it = indexNameMap.find(firstCol);
                if (it == indexNameMap.end()) {
                    
                    if (hasGroups) {  counts.push_back(thisGroupsCount);  }
                    indexNameMap[firstCol] = uniques;
                    totals.push_back(thisTotal);
                    total += thisTotal;
                    uniques++;
                    
                }else { error = true; m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + firstCol + ", sequence names must be unique. Please correct.\n"); }
            }
            in.close();

        }else if ((namefile == "") && (groupfile != "")) { //create count file from group only
            
            vector<string> names = groupMap->getNamesSeqs(); //only contains names from selectedGroups or all groups if selectedGroups is empty
            
            for (int i = 0; i < names.size(); i++) {
               if (m->getControl_pressed()) { break; }
                
                vector<countTableItem> abunds;
                string group = groupMap->getGroup(names[i]);
                int groupIndex = indexGroupMap[group];
                totalGroups[groupIndex]++;
                countTableItem thisAbund(1, groupIndex);
                abunds.push_back(thisAbund);
                
                map<string, int>::iterator it = indexNameMap.find(names[i]);
                if (it == indexNameMap.end()) {
                    
                    counts.push_back(abunds);
                    indexNameMap[names[i]] = uniques;
                    totals.push_back(1);
                    total++;
                    uniques++;
                    
                }else { error = true; m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + names[i] + ", sequence names must be unique. Please correct.\n"); }
            }
        }

        if (error) { m->setControl_pressed(true); }
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
int CountTable::readTable(string file, string format) {
    try {
        if (format == "fasta") {
            filename = file;
            ifstream in;
            util.openInputFile(filename, in);

            hasGroups = false;
            groups.clear();
            totalGroups.clear();
            indexGroupMap.clear();
            indexNameMap.clear();
            counts.clear();
            bool error = false;
            uniques = 0;
            total = 0;
            while (!in.eof()) {

                if (m->getControl_pressed()) { break; }

                Sequence seq(in); util.gobble(in);
                string name = seq.getName();
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + '\t' + toString(1) + "\n"); }

                map<string, int>::iterator it = indexNameMap.find(name);
                if (it == indexNameMap.end()) {
                    indexNameMap[name] = uniques;
                    totals.push_back(1);
                    total ++;
                    uniques++;
                }else {
                    error = true;
                    m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + name + ", sequence names must be unique. Please correct.\n");
                }
            }
            in.close();

            if (error) { m->setControl_pressed(true); }
        }else { m->mothurOut("[ERROR]: Unsupported format: " + format + ", please correct.\n"); m->setControl_pressed(true);  }

        return total;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "readTable");
        exit(1);
    }
}
/************************************************************/
int CountTable::readTable(string file, bool readGroups, bool mothurRunning) {
    try {
        return (readTable(file, readGroups, mothurRunning, nullVector));
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "readTable");
        exit(1);
    }
}
/************************************************************/
int CountTable::readTable(ifstream& in, bool readGroups, bool mothurRunning) {
    try {
        return (readTable(in, readGroups, mothurRunning, nullVector));
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "readTable");
        exit(1);
    }
}
/************************************************************/
bool CountTable::isCountTable(string file) {
    try {
        
        filename = file;
        ifstream in;
        util.openInputFile(filename, in);
        
        string headers = util.getline(in); util.gobble(in);
        
        if (headers[0] == '#') { //is this a count file in compressed form
            isCompressed = true;
            
            //read headers
            headers = util.getline(in); util.gobble(in); //gets compressed group name map line
            headers = util.getline(in); util.gobble(in);
        }
        vector<string> columnHeaders = util.splitWhiteSpace(headers);
        in.close();
        
        bool isCount = true;
        if (columnHeaders.size() >= 2) {
            vector<string> defaultHeaders = getHardCodedHeaders();
            if (defaultHeaders.size() >= 2) {
                if ((columnHeaders[0] != defaultHeaders[0]) && (columnHeaders[0] != "OTU_Label")) { isCount = false; }
                if (columnHeaders[1] != defaultHeaders[1]) { isCount = false; }
            }else { isCount = false; }
        }else { isCount = false; }
        
        return isCount;

    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "isCountTable");
        exit(1);
    }
}
/************************************************************/
int CountTable::readTable(string file, bool readGroups, bool mothurRunning, vector<string> selectedGroups) {
    try {
        filename = file;
        ifstream in;
        util.openInputFile(filename, in);
        
        readTable(in, readGroups, mothurRunning, selectedGroups);
        
        in.close();
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "readTable");
        exit(1);
    }
}
/************************************************************/
int CountTable::readTable(ifstream& in, bool readGroups, bool mothurRunning, vector<string> selectedGroups) {
    try {
        if (!readGroups) { selectedGroups.clear(); }

        string headers = util.getline(in); util.gobble(in);
        
        map<string, int> headerIndex2Group;
        //#1,F003D000	2,F003D002	3,F003D004	4,F003D006	5,F003D008	6,F003D142	7,F003D144	8,F003D146	9,F003D148	10,F003D150
        if (headers[0] == '#') { //is this a count file in compressed form
            isCompressed = true;
            
            //read headers
            headers = util.getline(in); util.gobble(in); //gets compressed group name map line
            headers = headers.substr(1);
            
            vector<string> groupNameHeaders = util.splitWhiteSpace(headers);
            
            for (int i = 0; i < groupNameHeaders.size(); i++) {
                string groupIndex = ""; string groupName = groupNameHeaders[i];
                util.splitAtComma(groupIndex, groupName);
                int a; util.mothurConvert(groupIndex, a);
                headerIndex2Group[groupName] = a-1;
            }
            
            headers = util.getline(in); util.gobble(in);
        }
        
        vector<string> columnHeaders = util.splitWhiteSpace(headers);
        
        int numGroupsInFile = 0;
        groups.clear();
        totalGroups.clear();
        indexGroupMap.clear();
        indexNameMap.clear();
        counts.clear();
        map<int, string> originalGroupIndexes;
        if ((columnHeaders.size() > 2) && readGroups) { hasGroups = true; numGroupsInFile = columnHeaders.size() - 2;  }
        
        set<string> setOfSelectedGroups;
        if (readGroups) {
            for (int i = 2; i < columnHeaders.size(); i++) {
                bool saveGroup = true;
                if (selectedGroups.size() != 0) {
                    if (!(util.inUsersGroups(columnHeaders[i], selectedGroups))) { saveGroup = false; }
                } //is this group in selected groups
                
                if (saveGroup) {
                    groups.push_back(columnHeaders[i]);
                    if (isCompressed) {
                        map<string, int>::iterator it = headerIndex2Group.find(columnHeaders[i]);
                        if (it != headerIndex2Group.end()) {
                            originalGroupIndexes[it->second] = columnHeaders[i];
                        }
                    }
                    else {  originalGroupIndexes[i-2] = columnHeaders[i];  }
                    totalGroups.push_back(0);
                    setOfSelectedGroups.insert(columnHeaders[i]);
                }
            }
        }
        
        //sort groups to keep consistent with how we store the groups in groupmap
        sort(groups.begin(), groups.end());
        for (int i = 0; i < groups.size(); i++) {  indexGroupMap[groups[i]] = i; }
        int numGroupsSelected = groups.size();

        bool error = false;
        string name;
        int thisTotal;
        uniques = 0;
        total = 0;
        while (!in.eof()) {

            if (m->getControl_pressed()) { break; }

            in >> name; util.gobble(in); in >> thisTotal; util.gobble(in);
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + '\t' + toString(thisTotal) + "\n"); }

            if ((thisTotal == 0) && !mothurRunning) { error=true; m->mothurOut("[ERROR]: Your count table contains a sequence named " + name + " with a total=0. Please correct.\n"); 
            }

            //if group info, then read it
            vector<int> groupCounts; groupCounts.resize(numGroupsSelected, 0);
            if (columnHeaders.size() > 2) { //file contains groups
                if (readGroups) { //user wants to save them
                    if (selectedGroups.size() != 0) {
                        //read this seqs groups abundances
                        thisTotal = 0;
                        if (isCompressed) {
                            string groupInfo = util.getline(in); util.gobble(in);
                            vector<string> groupNodes = util.splitWhiteSpace(groupInfo);
                            
                            vector<countTableItem> abunds;
                            for (int i = 0; i < groupNodes.size(); i++) { //for each non zero group count
                                string abund = groupNodes[i]; string thisgroup = "";
                                util.splitAtComma(thisgroup, abund);
                                int a; util.mothurConvert(abund, a);
                                int g; util.mothurConvert(thisgroup, g); g--;
                                string groupName = originalGroupIndexes[g]; //order of groups in file may not be sorted
                                
                                if (setOfSelectedGroups.count(groupName) != 0) { //we selected this group
                                    int thisIndex = indexGroupMap[groupName];
                                    countTableItem item(a, thisIndex);
                                    abunds.push_back(item);
                                    totalGroups[thisIndex] += a;
                                    thisTotal += a;
                                }
                            }
                            
                            groupCounts = expandAbunds(abunds);
                        }else {
                            for (int i = 0; i < numGroupsInFile; i++) {
                                int thisGroupAbund = 0;
                                in >> thisGroupAbund; util.gobble(in);
                                string groupName = originalGroupIndexes[i]; //order of groups in file may not be sorted
                                
                                if (setOfSelectedGroups.count(groupName) != 0) { //we selected this group
                                    int thisIndex = indexGroupMap[groupName];
                                    groupCounts[thisIndex] = thisGroupAbund;
                                    totalGroups[thisIndex] += thisGroupAbund;
                                    thisTotal += thisGroupAbund;
                                }
                            }
                        }
                    }else {
                        
                            if (isCompressed) {
                                string groupInfo = util.getline(in); util.gobble(in);
                                vector<string> groupNodes = util.splitWhiteSpace(groupInfo);
                                
                                vector<countTableItem> abunds;
                                for (int i = 0; i < groupNodes.size(); i++) { //for each non zero group count
                                    string abund = groupNodes[i]; string thisgroup = "";
                                    util.splitAtComma(thisgroup, abund);
                                    int a; util.mothurConvert(abund, a);
                                    int g; util.mothurConvert(thisgroup, g); g--;
                                    string groupName = originalGroupIndexes[g]; //order of groups in file may not be sorted
                                    int thisIndex = indexGroupMap[groupName];
                                    countTableItem item(a, thisIndex);
                                   
                                    abunds.push_back(item);
                                    totalGroups[thisIndex] += a;
                                }
                                
                                groupCounts = expandAbunds(abunds);
                            }
                            else {
                                for (int i = 0; i < numGroupsInFile; i++) {
                                    int thisIndex = indexGroupMap[originalGroupIndexes[i]];
                                    in >> groupCounts[thisIndex]; util.gobble(in);
                                    totalGroups[thisIndex] += groupCounts[thisIndex];
                                }
                            }
                        
                    }
                }else { //read and discard
                    util.getline(in); util.gobble(in);
                }
            }
            
            map<string, int>::iterator it = indexNameMap.find(name);
            if (it == indexNameMap.end()) {
                bool saveSeq = true;
                if (hasGroups && readGroups) {
                    vector<countTableItem> thisGroupsCount = compressAbunds(groupCounts);
                    if (thisGroupsCount.size() == 0) {  saveSeq = false; }
                    else { counts.push_back(thisGroupsCount); }
                }
                if (saveSeq) {
                    indexNameMap[name] = uniques;
                    totals.push_back(thisTotal);
                    total += thisTotal;
                    uniques++;
                }
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + name + ", sequence names must be unique. Please correct.\n");
            }
        }

        if (error) { m->setControl_pressed(true); }
        else { //check for zero groups
            if (hasGroups && readGroups) {
                for (int i = 0; i < totalGroups.size(); i++) {
                    if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n");
                        removeGroup(groups[i]);
                        i--;
                    }
                }
            }
        }
        
        //if the file has groups, but we didn't read them
        if (!readGroups) { hasGroups = false; }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "readTable");
		exit(1);
	}
}
/************************************************************/
int CountTable::readTable(string file, bool readGroups, bool mothurRunning, set<string> selectedSeqs) {
    try {
        filename = file;
        ifstream in;
        util.openInputFile(filename, in);
        
        string headers = util.getline(in); util.gobble(in);
        
        map<string, int> headerIndex2Group;
        //#1,F003D000	2,F003D002	3,F003D004	4,F003D006	5,F003D008	6,F003D142	7,F003D144	8,F003D146	9,F003D148	10,F003D150
        if (headers[0] == '#') { //is this a count file in compressed form
            isCompressed = true;
            
            //read headers
            headers = util.getline(in); util.gobble(in); //gets compressed group name map line
            headers = headers.substr(1);
            
            vector<string> groupNameHeaders = util.splitWhiteSpace(headers);
            
            for (int i = 0; i < groupNameHeaders.size(); i++) {
                string groupIndex = ""; string groupName = groupNameHeaders[i];
                util.splitAtComma(groupIndex, groupName);
                int a; util.mothurConvert(groupIndex, a);
                headerIndex2Group[groupName] = a-1;
            }
            
            headers = util.getline(in); util.gobble(in);
        }
        
        vector<string> columnHeaders = util.splitWhiteSpace(headers);
        
        int numGroupsInFile = 0;
        groups.clear();
        totalGroups.clear();
        indexGroupMap.clear();
        indexNameMap.clear();
        counts.clear();
        map<int, string> originalGroupIndexes;
        if ((columnHeaders.size() > 2) && readGroups) { hasGroups = true; numGroupsInFile = columnHeaders.size() - 2;  }
        
                
        if (readGroups) {
            for (int i = 2; i < columnHeaders.size(); i++) {
                groups.push_back(columnHeaders[i]);
                
                if (isCompressed) {
                    map<string, int>::iterator it = headerIndex2Group.find(columnHeaders[i]);
                    if (it != headerIndex2Group.end()) {
                        originalGroupIndexes[it->second] = columnHeaders[i];
                    }
                }else { originalGroupIndexes[i-2] = columnHeaders[i];  }
                totalGroups.push_back(0);
            }
        }
        
        //sort groups to keep consistent with how we store the groups in groupmap
        sort(groups.begin(), groups.end());
        for (int i = 0; i < groups.size(); i++) {  indexGroupMap[groups[i]] = i; }
        int numGroups = groups.size();
        
        bool error = false;
        string name;
        int thisTotal;
        uniques = 0;
        total = 0;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            in >> name; util.gobble(in); in >> thisTotal; util.gobble(in);
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + '\t' + toString(thisTotal) + "\n"); }
            
            if ((thisTotal == 0) && !mothurRunning) { error=true; m->mothurOut("[ERROR]: Your count table contains a sequence named " + name + " with a total=0. Please correct.\n");
            }
            
            vector<int> groupCounts; groupCounts.resize(numGroups, 0);
            if (columnHeaders.size() > 2) { //file contains groups
                if (readGroups) { //user wants to save them
                    if (isCompressed) {
                        string groupInfo = util.getline(in); util.gobble(in);
                        vector<string> groupNodes = util.splitWhiteSpace(groupInfo);
                        
                        vector<countTableItem> abunds;
                        for (int i = 0; i < groupNodes.size(); i++) { //for each non zero group count
                            string abund = groupNodes[i]; string thisgroup = "";
                            util.splitAtComma(thisgroup, abund);
                            int a; util.mothurConvert(abund, a);
                            int g; util.mothurConvert(thisgroup, g); g--;
                            string groupName = originalGroupIndexes[g]; //order of groups in file may not be sorted
                            int thisIndex = indexGroupMap[groupName];
                            countTableItem item(a, thisIndex);
                            
                            abunds.push_back(item);
                            totalGroups[thisIndex] += a;
                        }
                        
                        groupCounts = expandAbunds(abunds);
                    }
                    else {
                        for (int i = 0; i < numGroupsInFile; i++) { int thisIndex = indexGroupMap[originalGroupIndexes[i]]; in >> groupCounts[thisIndex]; util.gobble(in); totalGroups[thisIndex] += groupCounts[thisIndex]; }
                    }
                }else { util.getline(in); util.gobble(in); }//read and discard
            }
            
            map<string, int>::iterator it = indexNameMap.find(name);
            if (it == indexNameMap.end()) {
                bool saveSeq = true;
                if (selectedSeqs.count(name) == 0) { //don't save
                    saveSeq = false;
                }
                if (saveSeq) {
                    if (hasGroups && readGroups) {
                        vector<countTableItem> thisGroupsCount = compressAbunds(groupCounts);
                        counts.push_back(thisGroupsCount);
                    }
                    indexNameMap[name] = uniques;
                    totals.push_back(thisTotal);
                    total += thisTotal;
                    uniques++;
                }
            }else {
                error = true;
                m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + name + ", sequence names must be unique. Please correct.\n");
            }
        }
        in.close();
        
        if (error) { m->setControl_pressed(true); }
        else { //check for zero groups
            if (hasGroups && readGroups) {
                for (int i = 0; i < totalGroups.size(); i++) {
                    if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
                }
            }
        }
        
        //if the file has groups, but we didn't read them
         if (!readGroups) { hasGroups = false; }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "readTable");
        exit(1);
    }
}
/************************************************************/

int CountTable::zeroOutTable() {
  try {

		for(int i=0;i<counts.size();i++){
			for(int j=0;j<counts[0].size();j++){
                counts[j].clear();
			}
		}

		totals.assign(totals.size(), 0);

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "CountTable", "zeroOutTable");
		exit(1);
	}
}
/************************************************************/

int CountTable::clearTable() {
    try {
        hasGroups = false;
        total = 0;
        uniques = 0;
        groups.clear();
        counts.clear();
        totals.clear();
        totalGroups.clear();
        indexNameMap.clear();
        indexGroupMap.clear();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "clearTable");
        exit(1);
    }
}
/************************************************************/
//zeroed reads are not printed
vector<string> CountTable::printTable(string file) {
    try {
        
        //remove group if all reads are removed
        for (int i = 0; i < totalGroups.size(); i++) {
            if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
        }
        
        if (isCompressed) { return printCompressedTable(file); }
        
        ofstream out;
        util.openOutputFile(file, out);
        
        vector<string> namesInTable;
        
        if (total != 0) {
            printHeaders(out);
            
            map<int, string> reverse; //use this to preserve order
            for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
            
            for (int i = 0; i < totals.size(); i++) {
                if (totals[i] != 0) {

                    map<int, string>::iterator itR = reverse.find(i);
                
                    if (itR != reverse.end()) {
                        
                        namesInTable.push_back(itR->second);
                        
                        out << itR->second << '\t' << totals[i];
                        
                        if (hasGroups) { printGroupAbunds(out, i); }
                        
                        out << endl;
                    }
                }
            }
        }
        out.close();
        return namesInTable;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "printTable");
		exit(1);
	}
}
/************************************************************/
//zeroed reads are not printed
vector<string> CountTable::printNoGroupsTable(string file) {
    try {
        
        ofstream out;
        util.openOutputFile(file, out);
        
        vector<string> namesInTable;
        
        if (total != 0) {
            vector<string> headers = getHardCodedHeaders();
            out << headers[0] << '\t' << headers[1] << endl;
            
            map<int, string> reverse; //use this to preserve order
            for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
            
            for (int i = 0; i < totals.size(); i++) {
                if (totals[i] != 0) {

                    map<int, string>::iterator itR = reverse.find(i);
                
                    if (itR != reverse.end()) {
                        
                        namesInTable.push_back(itR->second);
                        
                        out << itR->second << '\t' << totals[i] << endl;
                    }
                }
            }
        }
        out.close();
        return namesInTable;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "printTable");
        exit(1);
    }
}
/************************************************************/
//zeroed reads are not printed
vector<string> CountTable::printTable(string file, bool compressedFormat) {
    try {
        
        //remove group if all reads are removed
        for (int i = 0; i < totalGroups.size(); i++) {
            if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
        }
        
        if (compressedFormat) { return printCompressedTable(file); }
        
        ofstream out;
        util.openOutputFile(file, out);
        
        vector<string> namesInTable;
        
        if (total != 0) {
            printHeaders(out);
            
            map<int, string> reverse; //use this to preserve order
            for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
            
            for (int i = 0; i < totals.size(); i++) {
                
                if (totals[i] != 0) {
                    
                    map<int, string>::iterator itR = reverse.find(i);
                    
                    if (itR != reverse.end()) {
                        namesInTable.push_back(itR->second);
                        
                        out << itR->second << '\t' << totals[i];
                        
                        if (hasGroups) { printGroupAbunds(out, i); }
                        
                        out << endl;
                    }
                }
            }
        }
        out.close();
        return namesInTable;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "printTable");
        exit(1);
    }
}
/************************************************************/
//zeroed seqs are not printed
vector<string> CountTable::printCompressedTable(string file, vector<string> groupsToPrint) {
    try {
        ofstream out;
        util.openOutputFile(file, out);
        
        vector<string> namesInTable;
        
        bool pickedGroups = false;
        set<int> selectedGroupsIndicies;
        if (groupsToPrint.size() != 0) { if (hasGroups) { pickedGroups = true; } } //if no groups selected, print all groups
        
        if (total != 0) {
            if (hasGroups) {
                
                map<int, string> reverse;
                for (map<string, int>::iterator it = indexGroupMap.begin(); it !=indexGroupMap.end(); it++) { reverse[it->second] = it->first; }
                
                map<int, string>::iterator it = reverse.begin();
                string group1Name = it->second;
                if (pickedGroups) { //find selected groups indicies
                    for (map<int, string>::iterator it = reverse.begin(); it != reverse.end(); it++) {
                        if (util.inUsersGroups(it->second, groupsToPrint)) { group1Name = it->second; break; }
                    }
                }
                
                out << "#Compressed Format: groupIndex,abundance. For example 1,6 would mean the read has an abundance of 6 for group " + group1Name + "." << endl;
                out << "#";
                
                for (map<int, string>::iterator it = reverse.begin(); it != reverse.end(); it++) {
                    if (pickedGroups) { //find selected groups indicies
                        if (util.inUsersGroups(it->second, groupsToPrint)) {
                            selectedGroupsIndicies.insert(it->first);
                            
                            out << it->first+1 << "," << it->second << "\t";
                        }
                    }else { out << it->first+1 << "," << it->second << "\t"; }
                }
                out << endl;
            }
            
            printHeaders(out, groupsToPrint);
            
            map<int, string> reverse; //use this to preserve order
            for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
            
            for (int i = 0; i < totals.size(); i++) {
                if (totals[i] != 0) {
                    if (pickedGroups) {
                        string groupOutput = "";
                        long long thisTotal = 0;
                        for (int j = 0; j < counts[i].size(); j++) {
                            if (selectedGroupsIndicies.count(counts[i][j].group) != 0) { //this is a group we want
                                groupOutput += '\t' + toString(counts[i][j].group+1) + ',' + toString(counts[i][j].abund);
                                thisTotal += counts[i][j].abund;
                            }
                        }
                        
                        if (thisTotal != 0) {
                            map<int, string>::iterator itR = reverse.find(i);
                            
                            if (itR != reverse.end()) {
                                namesInTable.push_back(itR->second);
                                
                                out << itR->second << '\t' << thisTotal << groupOutput << endl;
                            }
                        }
                    }
                    else {
                        map<int, string>::iterator itR = reverse.find(i);
                        
                        if (itR != reverse.end()) {
                            namesInTable.push_back(itR->second);
                            
                            out << itR->second << '\t' << totals[i];
                            if (hasGroups) {
                                for (int j = 0; j < counts[i].size(); j++) {
                                    out  << '\t' << counts[i][j].group+1 << ',' << counts[i][j].abund;
                                }
                            }
                            out << endl;
                        }
                    }
                }
            }
        }
        out.close();
        
        return namesInTable;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "printCompressedTable");
        exit(1);
    }
}
/************************************************************/
//returns index of countTableItem for group passed in. If group is not present in seq, returns index of next group or -1
int CountTable::find(int seq, int group, bool returnNext) {
    try {
        
        //if (!returnNext) { return find(seq, group); }
        int index = -1;
        
        for (int i = 0; i < counts[seq].size(); i++) {
            if (counts[seq][i].group >= group) { //found it or done looking
                
                if (counts[seq][i].group == group) { index = i;  }
                break;
            }
        }
        
        return index;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "find");
        exit(1);
    }
}/************************************************************/
//returns abundance of countTableItem for seq and group passed in. If group is not present in seq, returns 0
int CountTable::getAbund(int seq, int group) {
    try {
        int index = find(seq, group, false);
        
        if (index != -1) { //this seq has a non zero abundance for this group
            return counts[seq][index].abund;
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "getAbund");
        exit(1);
    }
}
/************************************************************/
vector<int> CountTable::expandAbunds(vector<countTableItem>& items) {
    try {
        vector<int> abunds; abunds.resize(groups.size(), 0); //prefill with 0's
        
        for (int i = 0; i < items.size(); i++) { //for each non zero entry
            abunds[items[i].group] = items[i].abund; //set abund for group
        }
        
        return abunds;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "expandAbunds");
        exit(1);
    }
}
/************************************************************/
vector<int> CountTable::expandAbunds(int index) {
    try {
        vector<int> abunds; abunds.resize(groups.size(), 0); //prefill with 0's
        
        
        for (int i = 0; i < counts[index].size(); i++) { //for each non zero entry
            abunds[counts[index][i].group] = counts[index][i].abund; //set abund for group
        }
        
        return abunds;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "expandAbunds");
        exit(1);
    }
}
/************************************************************/
//assumes same order as groups
vector<countTableItem> CountTable::compressAbunds(vector<int> abunds) {
    try {
        vector<countTableItem> row;
        
        for (int i = 0; i < abunds.size(); i++) {
            if (abunds[i] != 0) {
                countTableItem thisAbund(abunds[i], i);
                row.push_back(thisAbund);
            }
        }
        
        return row;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "compressAbunds");
        exit(1);
    }
}
/************************************************************/
void CountTable::printGroupAbunds(ofstream& out, int index) {
    try {
        
        vector<int> abunds = expandAbunds(index);
        
        for (int i = 0; i < abunds.size(); i++) { out << '\t' << abunds[i]; }
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "printGroupAbunds");
        exit(1);
    }
}
/************************************************************/
vector<string> CountTable::printSortedTable(string file) {
    try {
        //remove group if all reads are removed
        for (int i = 0; i < totalGroups.size(); i++) {
            if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
        }
        
        ofstream out;
        util.openOutputFile(file, out);
        printHeaders(out);
        
        vector<string> namesInTable;
        
        for (map<string, int>::iterator it = indexNameMap.begin(); it !=indexNameMap.end(); it++) {
            string seqName = it->first;
            int index = it->second;
            
            if (totals[index] != 0) {
                namesInTable.push_back(seqName);
                
                out << seqName << '\t' << totals[index];
                if (hasGroups) {
                    printGroupAbunds(out, index);
                }
                out << endl;
            }
        }
        out.close();
        
        return namesInTable;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "printSortedTable");
        exit(1);
    }
}

/************************************************************/
vector<string> CountTable::getHardCodedHeaders() {
    try {
        vector<string> headers; headers.push_back("Representative_Sequence"); headers.push_back("total");
        return headers;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "printHeaders");
        exit(1);
    }
}
/************************************************************/
int CountTable::printHeaders(ofstream& out, vector<string> selectedGroups) {
    try {
        //remove group if all reads are removed
        for (int i = 0; i < totalGroups.size(); i++) {
            if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
        }
        
        bool pickedGroups = false;
        if (selectedGroups.size() != 0) { pickedGroups = true; }
        
        out << "Representative_Sequence\ttotal";
        if (hasGroups) {
            for (int i = 0; i < groups.size(); i++) {
                if (pickedGroups) {
                    if (util.inUsersGroups(groups[i], selectedGroups)) {  out << '\t' << groups[i]; }
                }
                else { out << '\t' << groups[i]; }
            }
        }
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
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
        }else {
            if (totals[it->second] != 0) {
                out << it->first << '\t' << totals[it->second];
                
                if (hasGroups) { printGroupAbunds(out, it->second); }
                
                out << endl;
            }
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "printSeq");
		exit(1);
	}
}
/************************************************************/
SharedRAbundVectors* CountTable::getShared(vector<string> selected, map<string, string>& seqNameToOtuName) {
    try {
        
        if (selected.size() == 0) {}
        else { setNamesOfGroups(selected); }
        
        return getShared(seqNameToOtuName);
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "getShared");
        exit(1);
    }
}
/************************************************************/
SharedRAbundVectors* CountTable::getShared(map<string, string>& seqNameToOtuName) {
    try {
        SharedRAbundVectors* lookup = new SharedRAbundVectors();
        
        if (hasGroups) {
            for (int i = 0; i < groups.size(); i++) { //create blank rabunds for each group
                SharedRAbundVector* thisGroupsRabund = new SharedRAbundVector();
                thisGroupsRabund->setGroup(groups[i]);
                lookup->push_back(thisGroupsRabund);
            }
            
            //generate OTULabels
            vector<string> otuNames;
            util.getOTUNames(otuNames, counts.size(), "Otu");
            
            //create name map for seq -> otuName for use by other commands with associated files
            for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) { seqNameToOtuName[it->first] = otuNames[it->second]; }
            
            //add each "otu"
            for (int i = 0; i < counts.size(); i++) {
                vector<int> abunds = expandAbunds(i);
                lookup->push_back(abunds, otuNames[i]);
            }
            
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }
        

        return lookup;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "getShared");
        exit(1);
    }
}
/************************************************************/
//group counts for a seq
vector<int> CountTable::getGroupCounts(string seqName) {
    try {
        vector<countTableItem> temp = getItems(seqName);
        return (expandAbunds(temp)); 
        
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getGroupCounts");
		exit(1);
	}
}
/************************************************************/
//group counts for a seq
vector<countTableItem> CountTable::getItems(string seqName) {
    try {
        vector<countTableItem> temp;
        if (hasGroups) {
            map<string, int>::iterator it = indexNameMap.find(seqName);
            if (it == indexNameMap.end()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (util.inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                }
                m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                temp = counts[it->second];
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n"); m->setControl_pressed(true); }
        
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
                m->mothurOut("[ERROR]: group " + groupName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                return totalGroups[it->second];
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }

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
                m->mothurOut("[ERROR]: group " + groupName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                map<string, int>::iterator it2 = indexNameMap.find(seqName);
                if (it2 == indexNameMap.end()) {
                    //look for it in names of groups to see if the user accidently used the wrong file
                    if (util.inUsersGroups(seqName, groups)) {
                        m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                    }
                    m->mothurOut("[ERROR]: seq " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
                }else {
                    return expandAbunds(it2->second)[it->second];
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }

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
                m->mothurOut("[ERROR]: " + groupName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                map<string, int>::iterator it2 = indexNameMap.find(seqName);
                if (it2 == indexNameMap.end()) {
                    //look for it in names of groups to see if the user accidently used the wrong file
                    if (util.inUsersGroups(seqName, groups)) {
                        m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                    }
                    m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
                }else {
                    int indexOfGroup = find(it2->second, it->second, false);
                    int oldCount = 0;
                    
                    if (indexOfGroup == -1) { //create item for this group
                        countTableItem newItem(num, it->second);
                        counts[it2->second].push_back(newItem);
                        sortRow(it2->second);
                    }else { //update total for group
                        oldCount = counts[it2->second][indexOfGroup].abund;
                        counts[it2->second][indexOfGroup].abund = num;
                    }
                    
                    totalGroups[it->second] += (num - oldCount);
                    total += (num - oldCount);
                    totals[it2->second] += (num - oldCount);
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }

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
        bool sanity = util.inUsersGroups(groupName, groups);
        if (sanity) { m->mothurOut("[ERROR]: " + groupName + " is already in the count table, cannot add again.\n"); m->setControl_pressed(true);  return 0; }

        groups.push_back(groupName);
        if (!hasGroups) { counts.resize(uniques);  }
        
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

        hasGroups = true;

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
            map<string, int>::iterator it;
            for (it = indexNameMap.begin(); it !=indexNameMap.end(); it++) { reverse[it->second] = it->first;  }
            
            it = indexGroupMap.find(groupName);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + groupName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
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
                    
                    if (m->getControl_pressed()) { break; }
                    
                    int indexOfGroup = -1; bool found = false;
                    for (int j = 0; j < counts[i].size(); j++) {
                        if (counts[i][j].group >= indexOfGroupToRemove) { //found it or done looking
                            
                            indexOfGroup = j;
                            if (counts[i][j].group == indexOfGroupToRemove) {   found = true; }
                            break;
                        }
                    }
                    
                    if (found) { //you have an abundance for this group
                        int num = counts[i][indexOfGroup].abund;
                        counts[i].erase(counts[i].begin()+indexOfGroup);
                        totals[i] -= num;
                        total -= num;
                        
                        if (totals[i] == 0) { //your sequences are only from the group we want to remove, then remove you.
                            counts.erase(counts.begin()+i);
                            totals.erase(totals.begin()+i);
                            uniques--;
                            i--;
                            if (i == -1) { i = 0; }
                            indexOfGroup = counts[i].size(); //don't adjust the the group indexes because we removed the read
                        }else { newIndexNameMap[reverse[thisIndex]] = i; }
                    }else { //you don't have this group, nothing to remove
                        
                        if (indexOfGroup == -1) { indexOfGroup = counts[i].size(); }
                        newIndexNameMap[reverse[thisIndex]] = i;
                    }
                    
                    for (int j = indexOfGroup; j < counts[i].size(); j++) { counts[i][j].group -= 1; }
                    
                    thisIndex++;
                }
                indexNameMap = newIndexNameMap;
                
                if (groups.size() == 0) { hasGroups = false; }
            }
        }else { m->mothurOut("[ERROR]: your count table does not contain group information, can not remove group " + groupName + ".\n"); m->setControl_pressed(true); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "removeGroup");
        exit(1);
    }
}
/***********************************************************************/
int CountTable::removeGroup(int minSize){
    try {
        
        if (hasGroups) {
            for (int i = 0; i < totalGroups.size(); i++) {
                if (totalGroups[i] < minSize) { removeGroup(groups[i]); }
            }
        }else { m->mothurOut("[ERROR]: your count table does not contain group information, can not remove groups.\n"); m->setControl_pressed(true); }
        
        return groups.size();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "removeGroup - minSize");
        exit(1);
    }
}
/************************************************************/
//vector of groups for the seq
vector<string> CountTable::getGroups(string seqName) {
    try {
        vector<string> thisGroups;
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
        }else {
            if (hasGroups) {
                int index = it->second;
                for (int i = 0; i < counts[index].size(); i++) {
                    thisGroups.push_back(groups[counts[index][i].group]);
                }
            }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }
        }

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
                if (util.inUsersGroups(oldSeqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + oldSeqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                }
            }
            m->mothurOut("[ERROR]: " + oldSeqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
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
                if (util.inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                }
            }
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
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
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true); return -1;
        }else {
            int diff = totals[it->second] - abund;
            totals[it->second] = abund;
            total-=diff;
        }

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "setNumSeqs");
        exit(1);
    }
}
/************************************************************/
int CountTable::zeroOutSeq(string seqName) {
    try {

        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true); return -1;
        }else {
            int abund = totals[it->second];
            totals[it->second] = 0;
            total-=abund;
            
            if (hasGroups) {
                int seqIndexIntoCounts = it->second;
                for (int i = 0; i < counts[seqIndexIntoCounts].size(); i++) {
                    totalGroups[counts[seqIndexIntoCounts][i].group] -= counts[seqIndexIntoCounts][i].abund;
                }
                counts[seqIndexIntoCounts].clear();
            }
        }

        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "zeroOutSeq");
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
                if (util.inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n");
                }
            }
            m->mothurOut("[ERROR]: " + seqName + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
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
            if (hasGroups) {  m->mothurOut("[ERROR]: Your count table has groups and I have no group information for " + seqName + ".\n");  m->setControl_pressed(true);  }
            indexNameMap[seqName] = uniques;
            totals.push_back(1);
            total++;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct.\n");  m->setControl_pressed(true);
        }

        return 1;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "push_back");
		exit(1);
	}
}
/************************************************************/
//
bool CountTable::inTable(string seqName) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it != indexNameMap.end()) { return true; }
        return false;
        
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "inTable");
        exit(1);
    }
}

/************************************************************/
//remove sequence
int CountTable::remove(string seqName) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it != indexNameMap.end()) {
            int seqIndexIntoCounts = it->second;
            uniques--;
            if (hasGroups){ //remove this sequences counts from group totals
                for (int i = 0; i < counts[seqIndexIntoCounts].size(); i++) {
                    totalGroups[counts[seqIndexIntoCounts][i].group] -= counts[seqIndexIntoCounts][i].abund;
                }
            }
            
            //save for later in case removing a group means we need to remove a seq.
            map<int, string> reverse;
            for (map<string, int>::iterator it2 = indexNameMap.begin(); it2 !=indexNameMap.end(); it2++) { reverse[it2->second] = it2->first;  }
            
            int newIndex = 0;
            map<string, int> newIndexNameMap;
            for (int i = 0; i < counts.size(); i++) {
                if (i == seqIndexIntoCounts) { }//you are the seq we are trying to remove
                else {   newIndexNameMap[reverse[i]] = newIndex; newIndex++;  }
            }
            indexNameMap = newIndexNameMap;

            counts.erase(counts.begin()+seqIndexIntoCounts);
            int thisTotal = totals[seqIndexIntoCounts];
            totals.erase(totals.begin()+seqIndexIntoCounts);
            total -= thisTotal;
            
            //remove group if all reads are removed
            for (int i = 0; i < totalGroups.size(); i++) {
                if (totalGroups[i] == 0) { m->mothurOut("\nRemoving group: " + groups[i] + " because all sequences have been removed.\n"); removeGroup(groups[i]); i--; }
            }
            
        }else {
            if (hasGroupInfo()) {
                //look for it in names of groups to see if the user accidently used the wrong file
                if (util.inUsersGroups(seqName, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seqName + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                }
            }
            m->mothurOut("[ERROR]: Your count table contains does not include " + seqName + ", cannot remove.\n");  m->setControl_pressed(true);
        }

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "remove");
		exit(1);
	}
}
/************************************************************/
//add seqeunce without group info
int CountTable::push_back(string seqName, int thisTotal) {
    try {
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if (hasGroups) {  m->mothurOut("[ERROR]: Your count table has groups and I have no group information for " + seqName + ".\n"); m->setControl_pressed(true);  }
            indexNameMap[seqName] = uniques;
            totals.push_back(thisTotal);
            total+=thisTotal;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct.\n");  m->setControl_pressed(true);
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
int CountTable::push_back(string seqName, vector<int> groupCounts, bool ignoreDup=false) {
    try {
        int thisTotal = 0;
        map<string, int>::iterator it = indexNameMap.find(seqName);
        if (it == indexNameMap.end()) {
            if ((hasGroups) && (groupCounts.size() != getNumGroups())) {  m->mothurOut("[ERROR]: Your count table has a " + toString(getNumGroups()) + " groups and " + seqName + " has " + toString(groupCounts.size()) + ", please correct.\n");  m->setControl_pressed(true);  }
            
            for (int i = 0; i < getNumGroups(); i++) {   totalGroups[i] += groupCounts[i];  thisTotal += groupCounts[i]; }
            if (hasGroups) {  counts.push_back(compressAbunds(groupCounts));  }
            indexNameMap[seqName] = uniques;
            totals.push_back(thisTotal);
            total+= thisTotal;
            uniques++;
        }else {
            if (ignoreDup) {
                m->mothurOut("[WARNING]: Your count table contains more than 1 sequence named " + seqName + ".  Mothur requires sequence names to be unique. I will only add it once.\n"); 
            }else {  m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct.\n");  m->setControl_pressed(true);  }
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
            if ((hasGroups) && (groupCounts.size() != getNumGroups())) {  m->mothurOut("[ERROR]: Your count table has a " + toString(getNumGroups()) + " groups and " + seqName + " has " + toString(groupCounts.size()) + ", please correct.\n");  m->setControl_pressed(true);  }

            for (int i = 0; i < getNumGroups(); i++) {   totalGroups[i] += groupCounts[i];  thisTotal += groupCounts[i]; }
            if (hasGroups) {  counts.push_back(compressAbunds(groupCounts));  }
            indexNameMap[seqName] = uniques;
            totals.push_back(thisTotal);
            total+= thisTotal;
            uniques++;
        }else {
            m->mothurOut("[ERROR]: Your count table contains more than 1 sequence named " + seqName + ", sequence names must be unique. Please correct.\n");  m->setControl_pressed(true);
        }

        return thisTotal;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "push_back");
		exit(1);
	}
}
/************************************************************/
//returns size of smallest group. If no groups, returns total num seqs (includes non uniques)
int CountTable::getNumSeqsSmallestGroup() {
    try {
        int smallestGroupSize = MOTHURMAX;
        
        if (hasGroups) {
            for (int i = 0; i < totalGroups.size(); i++) {
                if (totalGroups[i] < smallestGroupSize) { smallestGroupSize = totalGroups[i]; }
            }
        }
        else { return total; }
        
        return smallestGroupSize;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "getNumSeqsSmallestGroup");
        exit(1);
    }
}

/************************************************************/
//create ListVector from uniques
ListVector CountTable::getListVector() {
    try {
        ListVector list(indexNameMap.size());
        for (map<string, int>::iterator it = indexNameMap.begin(); it != indexNameMap.end(); it++) {
            if (m->getControl_pressed()) { break; }
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
//returns the names of all unique sequences in file mapped to their seqCounts
map<string, int> CountTable::getNameMap(string group) {
    try {
        map<string, int> names;
        
        if (hasGroups) {
            map<string, int>::iterator it = indexGroupMap.find(group);
            if (it == indexGroupMap.end()) {
                m->mothurOut("[ERROR]: " + group + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                for (map<string, int>::iterator it2 = indexNameMap.begin(); it2 != indexNameMap.end(); it2++) {
                    int abund = getAbund(it2->second, it->second);
                    if (abund != 0) {  names[it2->first] = abund; }
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }

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
                m->mothurOut("[ERROR]: " + group + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                for (map<string, int>::iterator it2 = indexNameMap.begin(); it2 != indexNameMap.end(); it2++) {
                    if (getAbund(it2->second, it->second) != 0) {  names.push_back(it2->first); }
                }
            }
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }

        return names;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "getNamesOfSeqs");
		exit(1);
	}
}
/************************************************************/
//returns the names of all unique sequences in file
vector<string> CountTable::getNamesOfSeqs(vector<string> chosenGroups) {
    try {
        vector<string> names;
        if (hasGroups) {
            set<string> uniqueNames;
            for (int i = 0; i < chosenGroups.size(); i++) {
                vector<string> namesFromThisGroup = getNamesOfSeqs(chosenGroups[i]);
                for (int j = 0; j < namesFromThisGroup.size(); j++) { uniqueNames.insert(namesFromThisGroup[j]);  }
            }
            
            //only adds names once. seqs are likely present in more than one group, but we only want to enter them once
            for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) { names.push_back(*it); }
            
        }else{  m->mothurOut("[ERROR]: Your count table does not have group info. Please correct.\n");  m->setControl_pressed(true); }
        
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
                if (util.inUsersGroups(seq1, groups)) {
                    m->mothurOut("[WARNING]: Your group or design file contains a group named " + seq1 + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                }
            }
            m->mothurOut("[ERROR]: " + seq1 + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
        }else {
            map<string, int>::iterator it2 = indexNameMap.find(seq2);
            if (it2 == indexNameMap.end()) {
                if (hasGroupInfo()) {
                    //look for it in names of groups to see if the user accidently used the wrong file
                    if (util.inUsersGroups(seq2, groups)) {
                        m->mothurOut("[WARNING]: Your group or design file contains a group named " + seq2 + ".  Perhaps you are used a group file instead of a design file? A common cause of this is using a tree file that relates your groups (created by the tree.shared command) with a group file that assigns sequences to a group.\n"); 
                    }
                }
                m->mothurOut("[ERROR]: " + seq2 + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
            }else {
                if (hasGroupInfo()) { //if no group data then counts are empty
                    //merge data
                    vector<int> countsSeq1 = expandAbunds(it->second);
                    vector<int> countsSeq2 = expandAbunds(it2->second);
                
                    for (int i = 0; i < groups.size(); i++) { countsSeq1[i] += countsSeq2[i]; }
                
                    counts[it->second] = compressAbunds(countsSeq1);
                }
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
            push_back(names[i], thisCounts, false);
        }
        
        isCompressed = ct->isTableCompressed();

        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "CountTable", "copy");
		exit(1);
	}
}
/***********************************************************************/

int CountTable::sortCountTable(){
    try {
        
        //sorts each rows abunds by group
        //counts[i] = (1,4),(1,2),(3,7) -> (1,2),(1,4),(3,7)
        for (int i = 0; i < counts.size(); i++) {  sort(counts[i].begin(), counts[i].end(), compareGroups); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "sortCountTable");
        exit(1);
    }
}
/***********************************************************************/

int CountTable::sortRow(int index){
    try {
        
        //saves time in getSmallestCell, by making it so you dont search the repeats
        sort(counts[index].begin(), counts[index].end(), compareGroups);
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "CountTable", "sortRow");
        exit(1);
    }
}

/************************************************************/
