#ifndef CONCENSUSCOMMAND_H
#define CONCENSUSCOMMAND_H
/*
 *  consensuscommand.h
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
	ConcensusCommand(string);	
	~ConcensusCommand();
	int execute();
	void help();	
	
private:
	GlobalData* globaldata;
	SharedUtil* util;
	vector<Tree*> t;
	Tree* consensusTree;
	bool abort;
	vector<string> treeSet;		//set containing all members of the tree to start recursion.  filled in getSets().
	map< vector<string>, int > nodePairs;  //<map of possible combinations these combos are the pcounts or descendants info, to how many times they occured
										// ie. combos FI and EGK would create nodePairs[vector containing F and I] = 1; nodePairs[vector containing E, G and K] = 1
										// if you saw the combo FI again in another tree you would then update nodePairs[vector containing F and I] = 2;
										// requires vectors to be sorted to find key.
	map< vector<string>, int > nodePairsInTree;
	map<string, int>::iterator it;
	map< vector<string>, int>::iterator it2;
	string outputFile, notIncluded, filename;
	ofstream out, out2;
	int numNodes, numLeaves, count;  //count is the next available spot in the tree vector
									 	
	void getSets();
	vector<string> getNextAvailableSet(vector<string>);  //gets next largest and highest rated set that is a subset of the set passed in.
	vector<string> getRestSet(vector<string>, vector<string>);
	bool isSubset(vector<string>, vector<string>); 
	int findSpot(string); 
	int buildConcensusTree(vector<string>);
	
};

#endif

