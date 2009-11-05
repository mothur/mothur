/*
 *  doTaxonomy.h
 *  
 *
 *  Created by Pat Schloss on 6/17/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"

/**************************************************************************************************/

struct TaxNode {
	vector<string> accessions;
	map<string, int> children;
	int parent, childNumber, level;
	string name, heirarchyID;
	
	TaxNode(string n) : name(n), level(0), parent(-1) {		}
	TaxNode(){}
};

/**************************************************************************************************/

class PhyloTree {

public:
	PhyloTree();
	void addSeqToTree(string, string);
	void assignHeirarchyIDs(int);
	void print(ofstream&);
private:
	string getNextTaxon(string&);
	vector<TaxNode> tree;
	void print(int, ofstream&);
	int numNodes;
	int numSeqs;
};

/**************************************************************************************************/

