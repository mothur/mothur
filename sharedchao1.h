#ifndef SHAREDCHAO1_H
#define SHAREDCHAO1_H
/*
 *  sharedchao1.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the Sharedchao1 estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/


class SharedChao1 : public Calculator  {
	
	public: 
		SharedChao1() : Calculator("sharedchao", 3, true) {};
		EstOutput getValues(SAbundVector*) {return data;};
		EstOutput getValues(vector<SharedRAbundVector*>);
	private:
		IntNode* f1root;
		IntNode* f2root;
		vector<IntNode*> f1leaves;
		vector<IntNode*> f2leaves;
		int numLeaves;
		int numNodes;

		void initialTree(int);  //builds trees structure with n leaf nodes initialized to 0.
		void setCoef(IntNode*, int);
		void updateTree(vector<int>); //take vector containing the abundance info. for a bin and updates trees.
		void updateBranchf1(IntNode*, vector<int>, int);  //pointer, vector of abundance values, index into vector
		void updateBranchf2(IntNode*, vector<int>, int);  //pointer, vector of abundance values, index into vector
		
		//for debugging
		void printTree();
		void printBranch(IntNode*);
};

/***********************************************************************/

#endif
