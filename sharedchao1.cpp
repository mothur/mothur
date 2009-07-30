/*
 *  sharedchao1.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedchao1.h"

/***********************************************************************/
EstOutput SharedChao1::getValues(vector<SharedRAbundVector*> shared){
	try {
		data.resize(1,0);		
		vector<int> temp; 
		int numGroups = shared.size();
		float Chao = 0.0; float leftvalue, rightvalue;
				
		// IntNode is defined in mothur.h
		// The tree used here is a binary tree used to represent the f1+++, f+1++, f++1+, f+++1, f11++, f1+1+... 
		// combinations required to solve the chao estimator equation for any number of groups.  Conceptually, think
		// of each node as having a 1 and a + value, or for f2 values a 2 and a + value, and 2 pointers to intnodes, and 2 coeffient values.
		// The coeffient value is how many times you chose branch 1 to get to that fvalue.
		// If you choose left you are selecting the 1 or 2 value and right means the + value.  For instance, to find
		// the number of bins that have f1+1+ you would start at the root, go left, right, left, and select the rightvalue.
		// the coeffient is 2.  Note: we only set the coeffient in f2 values.
		
		//create and initialize trees to 0.
		initialTree(numGroups);	
		
		
		for (int i = 0; i < shared[0]->size(); i++) {
			//get bin values and calc shared 
			bool sharedByAll = true;
			temp.clear();
			for (int j = 0; j < numGroups; j++) {
				temp.push_back(shared[j]->getAbundance(i));
				if (temp[j] == 0) { sharedByAll = false; }
			}
			
			//they are shared
			if (sharedByAll == true) { 
				//find f1 and f2values
				updateTree(temp);
			}
		}

			
		//calculate chao1, (numleaves-1) because numleaves contains the ++ values.
		bool bias;
		for(int i=0;i<numLeaves;i++){
			if (f2leaves[i]->lvalue == 0 || f2leaves[i]->rvalue == 0) { bias = true;}// break;}
		}

		if(bias){
			for (int i = 0; i < numLeaves; i++) {
				
				leftvalue = (float)(f1leaves[i]->lvalue * (f1leaves[i]->lvalue - 1)) / (float)((pow(2, (float)f2leaves[i]->lcoef)) * (f2leaves[i]->lvalue + 1));
				if (i != (numLeaves-1)) {
					rightvalue = (float)(f1leaves[i]->rvalue * (f1leaves[i]->rvalue - 1)) / (float)((pow(2, (float)f2leaves[i]->rcoef)) * (f2leaves[i]->rvalue + 1));
				}else{
					//add in sobs
					rightvalue = (float)(f1leaves[i]->rvalue);
				}
				Chao += leftvalue + rightvalue;
			}
		}
		else{
			
			for (int i = 0; i < numLeaves; i++) {
				
				leftvalue = (float)(f1leaves[i]->lvalue * f1leaves[i]->lvalue) / (float)((pow(2, (float)f2leaves[i]->lcoef)) * f2leaves[i]->lvalue);
				if (i != (numLeaves-1)) {
					rightvalue = (float)(f1leaves[i]->rvalue * f1leaves[i]->rvalue) / (float)((pow(2, (float)f2leaves[i]->rcoef)) * f2leaves[i]->rvalue);
				}else{
					//add in sobs
					rightvalue = (float)(f1leaves[i]->rvalue);
				}
				Chao += leftvalue + rightvalue;
			}
		}
		
		for (int i = 0; i < numNodes; i++) {
			delete f1leaves[i];
			delete f2leaves[i];
		}
		
		
		data[0] = Chao;
		return data;
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "getValues");
		exit(1);
	}
}

/***********************************************************************/
//builds trees structure with n leaf nodes initialized to 0.
void SharedChao1::initialTree(int n) {  
	try {
		// (2^n) / 2. Divide by 2 because each leaf node contains 2 values. One for + and one for 1 or 2.
		numLeaves = pow(2, (float)n) / 2;
		numNodes = 2*numLeaves - 1;
		int countleft = 0;
		int countright = 1;
		
		f1leaves.resize(numNodes);
		f2leaves.resize(numNodes);
		
		//initialize leaf values
		for (int i = 0; i < numLeaves; i++) {
			f1leaves[i] = new IntNode;
			f1leaves[i]->lvalue = 0;
			f1leaves[i]->rvalue = 0;
			f1leaves[i]->left = NULL;
			f1leaves[i]->right = NULL;
			
			f2leaves[i] = new IntNode;
			f2leaves[i]->lvalue = 0;
			f2leaves[i]->rvalue = 0;
			f2leaves[i]->left = NULL;
			f2leaves[i]->right = NULL;
		}
		
		//set pointers to children
		for (int j = numLeaves; j < numNodes; j++) {
			f1leaves[j] = new IntNode;
			f1leaves[j]->left = f1leaves[countleft];
			f1leaves[j]->right = f1leaves[countright];
						
			f2leaves[j] = new IntNode;
			f2leaves[j]->left = f2leaves[countleft];
			f2leaves[j]->right =f2leaves[countright];
						
			countleft = countleft + 2;
			countright = countright + 2;
		}
		
		//point to root
		f1root = f1leaves[numNodes-1];
		
		//point to root
		f2root = f2leaves[numNodes-1];
		
		//set coeffients
		setCoef(f2root, 0);
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "initialTree");
		exit(1);
	}
}

/***********************************************************************/
//take vector containing the abundance info. for a bin and updates trees.
void SharedChao1::updateTree(vector<int> bin) { 
	try {
		updateBranchf1(f1root, bin, 0);  
		updateBranchf2(f2root, bin, 0); 
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "updateTree");
		exit(1);
	}
}

/***********************************************************************/
void SharedChao1::updateBranchf1(IntNode* node, vector<int> bin, int index) {
	try {
		//if you have more than one group
		if (index == (bin.size()-1)) {
			if (bin[index] == 1) { node->lvalue++; node->rvalue++; }
			else { node->rvalue++;  }
		}else {
			if (bin[index] == 1) {
				//follow path as if you are 1
				updateBranchf1(node->left, bin, index+1);
			}
			//follow path as if you are +
			updateBranchf1(node->right, bin, index+1);
		}
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "updateBranchf1");		
		exit(1);
	}
}

/***********************************************************************/
void SharedChao1::updateBranchf2(IntNode* node, vector<int> bin, int index) {
	try {
		//if you have more than one group
		if (index == (bin.size()-1)) {
			if (bin[index] == 2) { node->lvalue++; node->rvalue++; }
			else { node->rvalue++;  }
		}else {
			if (bin[index] == 2) {
				//follow path as if you are 1
				updateBranchf2(node->left, bin, index+1);
			}
			//follow path as if you are +
			updateBranchf2(node->right, bin, index+1);
		}
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "updateBranchf2");	
		exit(1);
	}
}

/***********************************************************************/
void SharedChao1::setCoef(IntNode* node, int coef) {
	try {
		if (node->left != NULL) {
			setCoef(node->left, coef+1);
			setCoef(node->right, coef);
		}else {
			node->lcoef = coef+1;
			node->rcoef = coef;
		}
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "setCoef");	
		exit(1);
	}
}

/***********************************************************************/
//for debugging purposes
void SharedChao1::printTree() {
	
	mothurOut("F1 leaves"); mothurOutEndLine();
	printBranch(f1root);
	
	mothurOut("F2 leaves"); mothurOutEndLine();
	printBranch(f2root);


}
/*****************************************************************/
void SharedChao1::printBranch(IntNode* node) {
	try {
		
		// you are not a leaf
		if (node->left != NULL) {
			printBranch(node->left);
			printBranch(node->right);
		}else { //you are a leaf
			mothurOut(toString(node->lvalue)); mothurOutEndLine();
			mothurOut(toString(node->rvalue)); mothurOutEndLine();
		}
		
	}
	catch(exception& e) {
		errorOut(e, "SharedChao1", "printBranch");	
		exit(1);
	}
}

/*****************************************************************/




