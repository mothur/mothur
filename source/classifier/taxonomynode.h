#ifndef TAXONOMYNODE
#define TAXONOMYNODE

/*
 *  taxonomynode.h
 *  
 *
 *  Created by Pat Schloss on 7/8/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

/**************************************************************************************************/

#include "mothurout.h"
/**************************************************************************************************/

class TaxonomyNode {
	
public:
	TaxonomyNode();
	TaxonomyNode(string, int);
	void setName(string);
	string getName();


	void setParent(int);
	int getParent();
	
	void makeChild(string, int);
	map<string, int> getChildren();
	int getChildIndex(string);
	int	getNumKids();
	int getNumSeqs();
	void setTotalSeqs(int);
	int getLevel();
	
private:
	int parent;
	map<string, int> children;
	int numChildren;
	int level;
	
protected:
    MothurOut* m;
	int numSeqs;
	int totalSeqs;
	string name;
};

/**************************************************************************************************/

#endif
