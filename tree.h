#ifndef TREE_H
#define TREE_H

/*
 *  tree.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

using namespace std;

#include <string>
#include <iostream>
#include <vector>
#include "treenode.h"
#include "globaldata.hpp"

/* This class represents the treefile. */



class Tree {
	public: 
		Tree();  
		~Tree() {};
		
		Tree generateRandomTree();
		void createNewickFile();
		int getIndex(string);
		void setIndex(string, int);
		int getNumNodes() { return numNodes; }
		int getNumLeaves(){	return numLeaves;}
		vector<Node> tree;		//the first n nodes are the leaves, where n is the number of sequences.
		
	private:
		GlobalData* globaldata;
		int findRoot();  //return index of root node
		void printBranch(int);  //recursively print out tree
		int numNodes, numLeaves;
		ofstream out;
		string filename;
};




#endif