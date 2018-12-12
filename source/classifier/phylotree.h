#ifndef DOTAXONOMY_H
#define DOTAXONOMY_H

/*
 * phylotree.h
 *  
 *
 *  Created by Pat Schloss on 6/17/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "mothurout.h"
#include "utils.hpp"
#include "currentfile.h"

/**************************************************************************************************/

struct TaxNode {
	vector<string> accessions;	//names of seqs in this branch of tree
	map<string, int> children;  //childs name to index in tree
	int parent, childNumber, level;
	string name, heirarchyID;
	
	TaxNode(string n) : name(n), level(0), parent(-1) {		}
	TaxNode(){}
};

/**************************************************************************************************/

class PhyloTree {
    
#ifdef UNIT_TEST
    friend class TestPhyloTree;
#endif


public:
	PhyloTree();
	PhyloTree(string);  //pass it a taxonomy file and it makes the tree
	PhyloTree(ifstream&, string);  //pass it a taxonomy file and it makes the train.tree
	~PhyloTree() {};
	int addSeqToTree(string, string);
	void assignHeirarchyIDs(int);
	void printTreeNodes(string); //used by bayesian to save time
	vector<int> getGenusNodes();
	vector<int> getGenusTotals();	
	void setUp(string);  //used to create file needed for summary file if you use () constructor and add seqs manually instead of passing taxonomyfile
		
	TaxNode get(int i);				
	TaxNode get(string seqName);
    vector<TaxNode> getNodes(int); //returns vector of nodes at given level
	string getName(int i);			
	int getGenusIndex(string seqName);	
	string getFullTaxonomy(string);	 //pass a sequence name return taxonomy
    vector<string> getSeqs(string);	 //returns names of sequences in given taxonomy
	
	int getMaxLevel()		{	return maxLevel;	}
	int getNumSeqs()		{	return numSeqs;		}
	int getNumNodes()		{	return (int)tree.size();	}
	
	bool ErrorCheck(vector<string>);
	
private:
	string getNextTaxon(string&, string);
	void print(ofstream&, vector<TaxNode>&); //used to create static reference taxonomy file
	void fillOutTree(int, vector<TaxNode>&); //used to create static reference taxonomy file
	void binUnclassified(string);
	
	vector<TaxNode> tree;
	vector<int> genusIndex; //holds the indexes in tree where the genus level taxonomies are stored
	vector<int> totals; //holds the numSeqs at each genus level taxonomy
	map<string, int> name2Taxonomy;  //maps name to index in tree
    map<string, int> name2GenusNodeIndex;
	set<int> uniqueTaxonomies;  //map of unique taxonomies
	map<int, int> leafNodes; //used to create static reference taxonomy file
	//void print(int, ofstream&);
	int numNodes;
	int numSeqs;
	int maxLevel;
	bool calcTotals;
	MothurOut* m;
    CurrentFile* current;
    Utils util;
};

/**************************************************************************************************/

#endif


