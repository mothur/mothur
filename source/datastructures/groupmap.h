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
#include "utils.hpp"

/* This class is a representation of the groupfile.  It is used by all the shared commands to determine what group a 
	certain sequence belongs to. */

class GroupMap {
public:
	GroupMap() { m = MothurOut::getInstance(); groupFileName = ""; }
	GroupMap(string);
	~GroupMap();
    
    int getCopy(GroupMap*);
    
    int readMap();
	int readMap(vector<string> groups); //selected groups read in. If groups.size() == 0, all groups are read
    int readMap(string, vector<string> groups); //filename, selected groups. selected groups read in. If groups.size() == 0, all groups are read
    int readMap(string);
	int readDesignMap();
    int readDesignMap(string);
    
	int getNumGroups();
	bool isValidGroup(string);  //return true if string is a valid group
	string getGroup(string);
    vector<string> getGroups(string); //returns groups represented by the seqs passed in. Think column two from a namefile row (seq1,seq2,seq3,seq4,seq5) -> (group1,group2). seqs1,seq3 are from group1, seq2,seq4,seq5 are from group2.
    vector<string> getGroups(vector<string>); //returns groups represented by the seqs passed in. Think column two from a namefile row (seq1,seq2,seq3,seq4,seq5) stored as a vector of names -> (group1,group2). seqs1,seq3 are from group1, seq2,seq4,seq5 are from group2.
    int getNumSeqs(string, string); //list of seq names, group. returns number of seqs from group passed represented by the seqs passed in. Think column two from a namefile row (seq1,seq2,seq3,seq4,seq5), group1 -> 2. seqs1,seq3 are from group1, seq2,seq4,seq5 are from group2.
    int getNumSeqs(vector<string>, string); //vector of seq names, group. returns number of seqs from group passed represented by the seqs passed in. Think column two from a namefile row (seq1,seq2,seq3,seq4,seq5), group1 -> 2. seqs1,seq3 are from group1, seq2,seq4,seq5 are from group2.
    
	void setGroup(string, string);
	vector<string> getNamesOfGroups() {
		sort(namesOfGroups.begin(), namesOfGroups.end());
		groupIndex.clear();
		for (int i = 0; i < namesOfGroups.size(); i++) { groupIndex[namesOfGroups[i]] = i; }
		return namesOfGroups;
	}
    
    void removeGroups(vector<string> groups);
    
    vector<string> getNamesSeqs();
    vector<string> getNamesSeqs(string); //get names of seqs belonging to group passed in
    vector<string> getNamesSeqs(vector<string>); //get names of seqs belonging to the set of groups passed in
	void setNamesOfGroups(vector<string> sn) { namesOfGroups = sn; }
	int getNumSeqs()  {  return (int)groupmap.size();  }
    int getNumSeqs(string); //return the number of seqs in a given group
    int getNumSeqsSmallestGroup(); //returns size of smallest group
	
    int renameSeq(string, string);
    int addSeq(string name, string group);
    
    int print(string);
    int print(ofstream&);
    int print(ofstream&, vector<string>); //print certain groups
    
    map<string, int> groupIndex;  //groupname, vectorIndex in namesOfGroups. - used by collectdisplays and libshuff commands.
    
private:
	vector<string> namesOfGroups;
	MothurOut* m;
	string groupFileName;
    int index;
	map<string, string>::iterator it;
	void setNamesOfGroups(string); 
	map<string, string> groupmap; //sequence name and groupname
	map<string, int> seqsPerGroup;  //maps groupname to number of seqs in that group
    Utils util;
};

#endif
