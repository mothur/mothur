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
//doesn't use MPI ifdefs since only pid 0 uses this class
class PhyloSummary {

public:
	PhyloSummary(string);
	PhyloSummary(string, string);
	~PhyloSummary() { if (groupmap != NULL)  {  delete groupmap;  }  }
	
	int summarize(string);  //pass it a taxonomy file and a group file and it makes the tree
	int addSeqToTree(string, string);
	int addSeqToTree(string, vector<string>);
	void print(ofstream&);
	int getMaxLevel() { return maxLevel; }
	
private:
	string getNextTaxon(string&);
	vector<rawTaxNode> tree;
	void print(int, ofstream&);
	void assignRank(int);
	void readTreeStruct(ifstream&);
	GroupMap* groupmap;
	bool ignore;
	
	int numNodes;
	int numSeqs;
	int maxLevel;
	MothurOut* m;
	
	void removeConfidences(string&);
};

/**************************************************************************************************/

#endif


