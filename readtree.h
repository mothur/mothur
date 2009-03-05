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

#include "globaldata.hpp"
#include "utilities.hpp"
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
		
	protected:
		GlobalData* globaldata;
		int numNodes, numLeaves;
		
};

/****************************************************************************/

class ReadNewickTree : public ReadTree {
	
public:
	ReadNewickTree(string file) : treeFile(file) { openInputFile(file, filehandle); } 
	~ReadNewickTree() {};
	void read();
	
private:
	Tree* T;
	int readNewickInt(istream&, int&, Tree*);
	void readTreeString();
	void nexusTranslation();
	ifstream filehandle;
	string treeFile;
	string holder;
};

/****************************************************************************/

#endif
