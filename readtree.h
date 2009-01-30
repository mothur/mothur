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

#include <string>
#include <iostream>
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
		
		virtual void read() {};
		int readSpecialChar(istream&, char, string);
		int readNodeChar(istream& f);
		float readBranchLength(istream& f);

				
		Tree* getTree()  { return T; }
		
		int numNodes, numLeaves;
		GlobalData* globaldata;
		
	protected:
		Tree* T;
};

/****************************************************************************/

class ReadNewickTree : public ReadTree {
	
public:
	ReadNewickTree(string file) : treeFile(file) { openInputFile(file, filehandle); }
	~ReadNewickTree() {};
	void read();
	
private:
	int readNewickInt(istream&, int&, Tree*);
	ifstream filehandle;
	string treeFile;
};

/****************************************************************************/

#endif
