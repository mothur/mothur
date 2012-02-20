#ifndef SUFFIXNODES_H
#define SUFFIXNODES_H

/*
 *  SuffixNodes.h
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	There are two types of nodes in a suffix tree as I have implemented it.  First, there are the internal nodes that
 *	have children, these are the SuffixBranch objects.  There are also the terminal nodes, which are the suffixBranches.
 *  I divided them into two groups to save on memory.  A SuffixTree object will be a vector of SuffixNodes; therefore,
 *	the values of parentNode, children nodes, and suffix nodes are stored as ints that correspond to indices in the 
 *	vector
 *
 */

#include "mothur.h"
#include "mothurout.h"

//********************************************************************************************************************

class SuffixNode {
	
public:
	SuffixNode(int, int, int);
	virtual ~SuffixNode() {}
	virtual void print(string, int)	= 0;
	virtual void setChildren(char, int);
	virtual int getNumChildren();
	virtual void eraseChild(char);
	virtual void setSuffixNode(int);
	virtual int getSuffixNode();
	virtual int getChild(char);
	int getParentNode();
	void setParentNode(int);
	int getStartCharPos();
	void setStartCharPos(int start);
	int getEndCharPos();
	
protected:
	int parentNode;
	int startCharPosition;
	int endCharPosition;
	MothurOut* m;
};

//********************************************************************************************************************

class SuffixLeaf : public SuffixNode {	//	most of the methods are already set in the parent class
	
public:
	SuffixLeaf(int, int, int);		//	we just need to define a constructor and
	~SuffixLeaf() {}
	void print(string, int);		//	print method
};

//********************************************************************************************************************

class SuffixBranch : public SuffixNode {
	
public:
	SuffixBranch(int, int, int);
	~SuffixBranch() {}
	void print(string, int);		//	need a special method for printing the node because there are children
	void eraseChild(char);			//	need a special method for erasing the children
	void setChildren(char, int);	//	need a special method for setting children
	void setSuffixNode(int);		//	need a special method for setting the suffix node
	int getSuffixNode();			//	need a special method for returning the suffix node
	int getChild(char);				//	need a special method for return children
	
private:
	vector<int> childNodes;			//	a suffix branch is unique because it has children and a suffixNode.  The 
	int suffixNode;					//	are stored in a vector for super-fast lookup.  If the alphabet were bigger, this
};									//	might not be practical.  Since we only have 5 possible letters, it makes sense

//********************************************************************************************************************

#endif
