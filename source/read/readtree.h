#ifndef READTREE_H
#define READTREE_H
/*
 *  readtree.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "tree.h"
#include "counttable.h"
#include "utils.hpp"

#define MAX_LINE		513
#define SKIPLINE(f,c)	{while((c=f.get())!=EOF && ((c) != '\n')){}}

class Tree;

/****************************************************************************/

class ReadTree {
	public:
		ReadTree(); 
		virtual ~ReadTree() {};
		
		virtual int read(CountTable*) = 0;
		int readSpecialChar(istream&, char, string);
		int readNodeChar(istream& f);
		float readBranchLength(istream& f);
	
		vector<Tree*> getTrees() { return Trees; }
		int AssembleTrees();
		
	protected:
		vector<Tree*> Trees;
		CountTable* ct;
		int numNodes, numLeaves;
		MothurOut* m;
        Utils util;
		
		
};

/****************************************************************************/

class ReadNewickTree : public ReadTree {
	
public:
    ReadNewickTree(string file, vector<string> T) : treeFile(file), Treenames(T) { Utils util; util.openInputFile(file, filehandle); readOk = 0; }
	~ReadNewickTree() {};
	int read(CountTable*);
	
private:
	Tree* T;
	int readNewickInt(istream&, int&, Tree*, CountTable*);
	int readTreeString(CountTable*);
	string nexusTranslation(CountTable*);
	ifstream filehandle;
	string treeFile;
	string holder;
	int readOk;  // readOk = 0 means success, readOk = 1 means errors.
    vector<string> Treenames;
	
};

/****************************************************************************/

#endif
