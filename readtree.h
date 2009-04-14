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

using namespace std;

#include "mothur.h"
#include "globaldata.hpp"
#include "tree.h"

#define MAX_LINE		513
#define SKIPLINE(f,c)	{while((c=f.get())!=EOF && ((c) != '\n')){}}

class Tree;

/****************************************************************************/

class ReadTree {
	public:
		ReadTree(); 
		~ReadTree() {};
		
		virtual int read() = 0;
		int readSpecialChar(istream&, char, string);
		int readNodeChar(istream& f);
		float readBranchLength(istream& f);
		
	protected:
		GlobalData* globaldata;
		int numNodes, numLeaves;
		
		
};

/****************************************************************************/

class ReadNewickTree : public ReadTree {
	
public:
	ReadNewickTree(string file) : treeFile(file) { openInputFile(file, filehandle); readOk = 0; } 
	~ReadNewickTree() {};
	int read();
	
private:
	Tree* T;
	int readNewickInt(istream&, int&, Tree*);
	int readTreeString();
	void nexusTranslation();
	ifstream filehandle;
	string treeFile;
	string holder;
	int readOk;  // readOk = 0 means success, readOk = 1 means errors.
	
};

/****************************************************************************/

#endif
