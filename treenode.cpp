/*
 *  treenode.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "treenode.h"


/****************************************************************/	
Node::Node() {
	//initialize node
	name = "";
	branchLength = -1;
	parent = -1;
	lchild = -1;
	rchild = -1;
	length2leaf = 0.0;
	
}
/****************************************************************/
void Node::setName(string Name) {  name = Name; }
/****************************************************************/
void Node::setGroup(string groups)  { group =groups; }
/****************************************************************/
void Node::setBranchLength(float l) { branchLength = l; }
/****************************************************************/
void Node::setLengthToLeaves(float l) { length2leaf = l; }
/****************************************************************/
void Node::setParent(int p)  { parent = p; }
/****************************************************************/
void Node::setIndex(int i)  { vectorIndex = i; }
/****************************************************************/
void Node::setChildren(int lc, int rc) { lchild = lc; rchild = rc; }	//leftchild, rightchild
/****************************************************************/
string Node::getName() { return name; }
/****************************************************************/
string Node::getGroup() { return group; }
/****************************************************************/
float Node::getBranchLength() { return branchLength; }
/****************************************************************/
float Node::getLengthToLeaves() { return length2leaf; }
/****************************************************************/
int Node::getParent() { return parent; }
/****************************************************************/
int Node::getLChild() { return lchild; }
/****************************************************************/
int Node::getRChild() { return rchild; }
/****************************************************************/
int Node::getIndex() { return vectorIndex; }
/****************************************************************/
//to be used by printTree in the Tree class to print the leaf info			
void Node::printNode() {
	try{
		cout << parent << ' ' << lchild << ' ' << rchild << ' ' << group;
		//there is a branch length
		if (branchLength != -1) { 
			cout << ' ' << setprecision(4) << branchLength; 
		}
		cout << " |";
		map<string, int>::iterator it;
		for(it=pGroups.begin();it!=pGroups.end();it++){
			cout << ' ' << it->first << ':' << it->second;
		}
		cout << " |";
		for(it=pcount.begin();it!=pcount.end();it++){
			cout << ' ' << it->first << ':' << it->second;
		}
		cout << endl; 
		
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Node class Function printNode. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Node class function printNode. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/****************************************************************/	
