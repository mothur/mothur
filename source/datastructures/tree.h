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

#include "treenode.h"
#include "counttable.h"
#include "currentfile.h"
/* This class represents the treefile. */

class Tree {
public: 
	Tree(int, CountTable*, vector<string>&);
	Tree(CountTable*, vector<string>&);		
    Tree(CountTable*, vector< vector<double> >&, vector<string>&); //create tree from sim matrix
	~Tree();
	
    CountTable* getCountTable() { return ct; }
    vector<string> getTreeNames() { return Treenames; }
	void getCopy(Tree*);  //makes tree a copy of the one passed in.
    void getCopy(Tree* copy, bool); //makes a copy of the tree structure passed in, (just parents, children and br). Used with the Tree(TreeMap*) constructor. Assumes the tmap already has set seqs groups you want.  Used by subsample to reassign seqs you don't want included to group "doNotIncludeMe".
	void getSubTree(Tree*, vector<string>);  //makes tree a that contains only the names passed in.
    
    //this function takes the leaf info and populates the non leaf nodes
    int assembleTree();
    void assembleRandomTree(Utils*); //pass tree indexes in random order
    void assembleRandomUnifracTree(vector<int>); //pass nodes to swap in random order
    
	void createNewickFile(string);
	int getIndex(string);
	void setIndex(string, int);
	int getNumNodes() { return numNodes; }
	int getNumLeaves(){	return numLeaves; }
	map<string, int> mergeUserGroups(int, vector<string>);  //returns a map with a groupname and the number of times that group was seen in the children
	void printTree();
	void print(ostream&);
	void print(ostream&, string);
    void print(ostream&, map<string, string>);
	int findRoot();  //return index of root node
    
	vector<Node> tree;		//the first n nodes are the leaves, where n is the number of sequences.
	map< string, vector<int> > groupNodeInfo;	//maps group to indexes of leaf nodes with that group, different groups may contain same node because of names file.
    vector<int> getNodes(vector<string>); //return tree indexes of nodes for groups passed in
			
private:
	CountTable* ct;
	int numNodes, numLeaves;
	ofstream out;
	string filename;
	
    //map<string, string> names;
	map<string, int>::iterator it, it2;
	map<string, int> mergeGroups(int);  //returns a map with a groupname and the number of times that group was seen in the children
	map<string,int> mergeGcounts(int);
    map<string, int> indexes; //maps seqName -> index in tree vector
	
    int randomLabels(vector<int>& nodesToSwap);
    int swapLabels(int first, int second);
	void addNamesToCounts(map<string, string>);
	void randomTopology(Utils*);
	void randomLabels(vector<string>);
	void printBranch(int, ostream&, map<string, string>);  //recursively print out tree
    void printBranch(int, ostream&, string);
	int populateNewTree(vector<Node>&, int, int&);
	void printBranch(int, ostream&, string, vector<Node>&);
    void pruneNewTree(Tree* copy, vector<string> namesToInclude);
    
	MothurOut* m;
    vector<string> Treenames;
    Utils util;
};

#endif
