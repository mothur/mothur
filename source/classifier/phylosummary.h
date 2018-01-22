#ifndef RAWTRAININGDATAMAKER_H
#define RAWTRAININGDATAMAKER_H

/*
 *  rawTrainingDataMaker.h
 *  Mothur
 *
 *  Created by westcott on 4/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"
#include "groupmap.h"
#include "counttable.h"
#include "utils.hpp"

/**************************************************************************************************/

struct rawTaxNode {
	map<string, int> children;  //childs name to index in tree
	int parent, level;
	string name, rank;
	map<string, int> groupCount;
	int total;
	
	rawTaxNode(string n) : name(n), level(0), parent(-1), total(0) {}
	rawTaxNode(){}
};

/**************************************************************************************************/
class PhyloSummary {

public:
	PhyloSummary(GroupMap*, bool, int);
	PhyloSummary(string, GroupMap*, bool, int);
    PhyloSummary(CountTable*, bool, int);
	PhyloSummary(string, CountTable*, bool, int);
	~PhyloSummary() {}
	
	int summarize(string);  //pass it a taxonomy file and a group file and it makes the tree
	int addSeqToTree(string, string);
	int addSeqToTree(string, map<string, bool>);
	void print(ofstream&, string);
    void print(ofstream&, bool);
	int getMaxLevel() { return maxLevel; }
	
private:
	string getNextTaxon(string&);
	vector<rawTaxNode> tree;
	void print(int, ofstream&, string);
    void print(int, ofstream&, bool);
	void assignRank(int);
    string getTaxons(vector<int> indexes, int index, int i, string&);
	void readTreeStruct(ifstream&);
    string findTaxon(string);
	GroupMap* groupmap;
    CountTable* ct;
	bool ignore, relabund;
	
	int numNodes, printlevel;
	int numSeqs;
	int maxLevel;
	MothurOut* m;
    Utils util;
};

/**************************************************************************************************/

#endif


