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

#define MAX_LINE		513
#define SKIPLINE(f,c)	{while((c=f.get())!=EOF && ((c) != '\n')){}}

class Tree;

/****************************************************************************/

class ReadTree {
	public:
		ReadTree(); 
		virtual ~ReadTree() {};
		
		virtual int read(TreeMap*) = 0;
		int readSpecialChar(istream&, char, string);
		int readNodeChar(istream& f);
		float readBranchLength(istream& f);
	
		vector<Tree*> getTrees() { return Trees; }
		int AssembleTrees(map<string, string>);
		
	protected:
		vector<Tree*> Trees;
		TreeMap* treeMap;
		int numNodes, numLeaves;
		MothurOut* m;
		
		
};

/****************************************************************************/

class ReadNewickTree : public ReadTree {
	
public:
	ReadNewickTree(string file) : treeFile(file) { m->openInputFile(file, filehandle); readOk = 0; } 
	~ReadNewickTree() {};
	int read(TreeMap*);
	
private:
	Tree* T;
	int readNewickInt(istream&, int&, Tree*, TreeMap*);
	int readTreeString(TreeMap*);
	string nexusTranslation(TreeMap*);
	ifstream filehandle;
	string treeFile;
	string holder;
	int readOk;  // readOk = 0 means success, readOk = 1 means errors.
	
};

/****************************************************************************/

#endif
