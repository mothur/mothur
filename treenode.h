#ifndef TREENODE_H
#define TREENODE_H

/*
 *  treenode.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

using namespace std;

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

/* This class represents a node on a tree. */


class Node  {
	public:
		Node();  //pass it the sequence name
		~Node() {};
		
		void setName(string);
		void setGroup(string);  //non leaf nodes will belong to multiple groups, leaf nodes will only belong to one.
		void setBranchLength(float);
		void setParent(int);
		void setChildren(int, int);		//leftchild, rightchild
		void setIndex(int);
		
		string getName();
		vector<string> getGroup();    //leaf nodes will only have 1 group, but branch nodes may have multiple groups.
		float getBranchLength();
		int getParent();
		int getLChild();
		int getRChild();
		int getIndex();
		
		void printNode(ostream&);   //prints out the name and the branch length
		
	private:
		string			name;
		vector<string>	group;
		float			branchLength;
		int				parent;
		int				lchild;
		int				rchild;
		int				vectorIndex;
};		

#endif
