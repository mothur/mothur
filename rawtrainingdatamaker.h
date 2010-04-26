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

/**************************************************************************************************/

struct rawTaxNode {
	map<string, int> children;  //childs name to index in tree
	int parent, level;
	string name, rank;
	
	rawTaxNode(string n) : name(n), level(0), parent(-1) {		}
	rawTaxNode(){}
};

/**************************************************************************************************/

class RawTrainingDataMaker {

public:
	RawTrainingDataMaker();
	RawTrainingDataMaker(string);  //pass it a taxonomy file and it makes the tree
	~RawTrainingDataMaker() {};
	int addSeqToTree(string, string);
	void assignRank(int);
	void print(ofstream&);
	
private:
	string getNextTaxon(string&);
	vector<rawTaxNode> tree;
	void print(int, ofstream&);
	int numNodes;
	int numSeqs;
	int maxLevel;
	//map<string, string> sanityCheck;
	MothurOut* m;
};

/**************************************************************************************************/

#endif


