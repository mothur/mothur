#ifndef TREE_H
#define TREE_H

/*
 *  tree.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

using namespace std;

#include <string>
#include <iostream>
#include <vector>

struct Node  {
		string	name;
		string	group;
		float	branchLength;
		Node*	parent;
		Node*	lchild;
		Node*	rchild;
};		



class Tree {
	public: 
		Tree();
		~Tree();
		
		Node* getParent(Node);
		Node* getLChild(Node);
		Node* getRChild(Node);
		
		void setParent(Node);
		void setLChild(Node);
		void setRChild(Node);
		
		
		Tree generateRandomTree();
		
		vector<Node> leaves;		//gives you easy access to the leaves of the tree to generate the parsimony score
		
	private:
};




#endif