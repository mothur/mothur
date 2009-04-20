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

#include "treenode.h"
#include "globaldata.hpp"

/* This class represents the treefile. */

class Tree {
public: 
	Tree();		//to generate a tree from a file
	~Tree();
	
	
	void getCopy(Tree*);  //makes tree a copy of the one passed in.
	void assembleRandomTree();
	void assembleRandomUnifracTree(vector<string>);
	void assembleRandomUnifracTree(string, string);
	void createNewickFile(string);
	int getIndex(string);
	void setIndex(string, int);
	int getNumNodes() { return numNodes; }
	int getNumLeaves(){	return numLeaves; }
	map<string, int> mergeUserGroups(int, vector<string>);  //returns a map with a groupname and the number of times that group was seen in the children
	void printTree();
	void print(ostream&);
	
	//this function takes the leaf info and populates the non leaf nodes
	void assembleTree();		
	
	vector<Node> tree;		//the first n nodes are the leaves, where n is the number of sequences.
private:
	GlobalData* globaldata;
	int numNodes, numLeaves;
	ofstream out;
	string filename;
	
	map<string, int>::iterator it, it2;
	map<string, int> mergeGroups(int);  //returns a map with a groupname and the number of times that group was seen in the children
	
	map<string,int> mergeGcounts(int);
	void randomTopology();
	void randomBlengths();
	void randomLabels(vector<string>);
	void randomLabels(string, string);
	int findRoot();  //return index of root node
	void printBranch(int, ostream&);  //recursively print out tree
};

#endif
