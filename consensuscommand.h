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

class ConcensusCommand : public Command {
	
public:
	ConcensusCommand(string);	
	ConcensusCommand();
	~ConcensusCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "concensus";	}
	string getCommandCategory()		{ return "Hidden";		}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	SharedUtil* util;
	vector<Tree*> t;
	Tree* consensusTree;
	bool abort;
	vector<string> treeSet;		//set containing all members of the tree to start recursion.  filled in getSets().
	map< vector<string>, int > nodePairs;  //<map of possible combinations these combos are the pcounts or descendants info, to how many times they occured
										// ie. combos FI and EGK would create nodePairs[vector containing F and I] = 1; nodePairs[vector containing E, G and K] = 1
										// if you saw the combo FI again in another tree you would then update nodePairs[vector containing F and I] = 2;
										// requires vectors to be sorted to find key.
	map< vector<string>, vector< vector<string> > > bestSplit;  //maps a group to its best split
	map< vector<string>, int > nodePairsInitialRate;
	map< vector<string>, int > nodePairsInTree;
	map<string, int>::iterator it;
	map< vector<string>, int>::iterator it2;
	string outputFile, notIncluded, filename;
	ofstream out, out2;
	int numNodes, numLeaves, count;  //count is the next available spot in the tree vector
	vector<string> outputNames;
									 	
	int getSets();
	int getSubgroupRating(vector<string>);
	vector<string> getSmallest(map< vector<string>, int>);
	vector<string> getNextAvailableSet(vector<string>, vector<string>&);  
	vector<string> getRestSet(vector<string>, vector<string>);
	bool isSubset(vector<string>, vector<string>); 
	int findSpot(string); 
	int buildConcensusTree(vector<string>);
	
};

#endif

