#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

/*
 *  suffixtree.h
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is my half-assed attempt to implement a suffix tree.  This is a cobbled together algorithm using materials that
 *	I found at http://marknelson.us/1996/08/01/suffix-trees/ and:
 *
 *		Ukkonen E. (1995). On-line construction of suffix trees. Algorithmica 14 (3): 249--260
 *		Gusfield, Dan (1999). Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology. 
 *			USA: Cambridge University Press
 *
 *	The Ukkonen paper is the seminal paper describing the on-line method of constructing a suffix tree.
 *
 *	I have chosen to store the nodes of the tree as a vector of pointers to SuffixNode objects.  The root is stored at
 *	nodeVector[0].  Each tree also stores the sequence name and the string that corresponds to the actual sequence. 
 *	Finally, this class provides a way of counting the number of suffixes that are needed in one tree to generate a new
 *	sequence (countSuffixes).  This method is used to determine similarity between sequences and was inspired by the
 *	article and Perl source code provided at http://www.ddj.com/web-development/184416093.
 *
 */

#include "mothur.h"

class SuffixNode;

//********************************************************************************************************************

class SuffixTree {
	
public:
	SuffixTree();
	~SuffixTree();

	void loadSequence(Sequence);
	string getSeqName();
	void print();	
	int countSuffixes(string, int&);
	int countSuffixes(string);	

private:	
	void addPrefix(int);
	void canonize();
	int split(int, int);
	void makeSuffixLink(int&, int);
	bool isExplicit(){	return activeStartPosition > activeEndPosition;	}
	
	int activeStartPosition;
	int activeEndPosition;
	
	vector<SuffixNode*> nodeVector;
	int root;
	int activeNode;
	int nodeCounter;
	string seqName;
	string sequence;
	MothurOut* m;
	
};

//********************************************************************************************************************

#endif
