#ifndef GROUPMAP_H
#define GROUPMAP_H
/*
 *  groupmap.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 12/1/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"

/* This class is a representation of the groupfile.  It is used by all the shared commands to determine what group a 
	certain sequence belongs to. */

class GroupMap {
public:
	GroupMap() { m = MothurOut::getInstance(); }
	GroupMap(string);
	~GroupMap();
	int readMap();
    int readMap(string);
	int readDesignMap();
    int readDesignMap(string);
	int getNumGroups();
	bool isValidGroup(string);  //return true if string is a valid group
	string getGroup(string);
	void setGroup(string, string);
	vector<string> getNamesOfGroups() {
		sort(namesOfGroups.begin(), namesOfGroups.end());
		groupIndex.clear();
		for (int i = 0; i < namesOfGroups.size(); i++) { groupIndex[namesOfGroups[i]] = i; }
		return namesOfGroups;
	}
    vector<string> getNamesSeqs();
	void setNamesOfGroups(vector<string> sn) { namesOfGroups = sn; }
	int getNumSeqs()  {  return groupmap.size();  }
	vector<string> getNamesSeqs(vector<string>); //get names of seqs belonging to a group or set of groups
	int getNumSeqs(string); //return the number of seqs in a given group
    int getCopy(GroupMap*);
    
    
    map<string, int> groupIndex;  //groupname, vectorIndex in namesOfGroups. - used by collectdisplays and libshuff commands.
    
private:
	vector<string> namesOfGroups;
	MothurOut* m;
	ifstream fileHandle;
	string groupFileName;
    int index;
	map<string, string>::iterator it;
	void setNamesOfGroups(string); 
	map<string, string> groupmap; //sequence name and groupname
	map<string, int> seqsPerGroup;  //maps groupname to number of seqs in that group
};

#endif
