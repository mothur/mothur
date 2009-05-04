#ifndef CONCENSUSCOMMAND_H
#define CONCENSUSCOMMAND_H
/*
 *  concensuscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/29/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "tree.h"
#include "treemap.h"
#include "sharedutilities.h"
	
class GlobalData;

class ConcensusCommand : public Command {
	
public:
	ConcensusCommand();	
	~ConcensusCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	SharedUtil* util;
	vector<Tree*> t;
	Tree* concensusTree;
	vector< map<string, int> > nodePairs;  //<maps a pair of nodes joined, number of times that pair occurred>	-one entry in vector for each internal node.
										// i.e. if node 7's child pairs are 1,2 ten times and 1,3 20 times then the map would be [12, 10] and [13, 20];
	map<string, int>::iterator it;
	map<string, int>::iterator it2;
	string outputFile, notIncluded;
	ofstream out, out2;
	int numNodes, numLeaves;
	
	void getNames(string, string&, string&);
};

#endif

