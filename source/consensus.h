#ifndef CONCENSUS_H
#define CONCENSUS_H
/*
 *  consensus.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/29/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */


#include "tree.h"
#include "treemap.h"

//NOTE: This class assumes all leaf nodes have 1 member.  
//      Mothur does allow for names files with trees which would make a tree with multiple members at one leaf.  
//      This class is currently only called internally by commands that have leaf node containing only 1 member.
//      But if in the future, this changes things will need to be reworked in getSets and buildConsensus.


class Consensus {
	
public:
	Consensus() { m = MothurOut::getInstance(); }	
	~Consensus() {}
	
    Tree* getTree(vector<Tree*>&);
		
private:
    MothurOut* m;
    Tree* consensusTree;
    
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
	int numNodes, numLeaves, count, numTrees;  //count is the next available spot in the tree vector
	vector<string> outputNames;
									 	
	int getSets(vector<Tree*>&);
	int getSubgroupRating(vector<string>);
	vector<string> getSmallest(map< vector<string>, int>);
	vector<string> getNextAvailableSet(vector<string>, vector<string>&);  
	vector<string> getRestSet(vector<string>, vector<string>);
	bool isSubset(vector<string>, vector<string>); 
	int findSpot(string); 
	int buildConsensusTree(vector<string>);
    int printSetsInfo();
	
};

#endif

