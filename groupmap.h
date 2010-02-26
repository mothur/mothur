#ifndef GROUPMAP_H
#define GROUPMAP_H
/*
 *  groupmap.h
 *  Dotur
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
	GroupMap() {};
	GroupMap(string);
	~GroupMap();
	int readMap();
	int getNumGroups();
	bool isValidGroup(string);  //return true if string is a valid group
	string getGroup(string);
	void setGroup(string, string);
	vector<string> namesOfGroups;
	map<string, int> groupIndex;  //groupname, vectorIndex in namesOfGroups. - used by collectdisplays and libshuff commands.
	int getNumSeqs()  {  return groupmap.size();  }
	vector<string> getNamesSeqs();
	int getNumSeqs(string); //return the number of seqs in a given group
			
private:
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
