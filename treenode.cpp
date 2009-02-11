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
	
}
/****************************************************************/
void Node::setName(string Name) {  name = Name; }
/****************************************************************/
void Node::setGroup(string groups)  { group =groups; }
/****************************************************************/
void Node::setBranchLength(float l) { branchLength = l; }
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
int Node::getParent() { return parent; }
/****************************************************************/
int Node::getLChild() { return lchild; }
/****************************************************************/
int Node::getRChild() { return rchild; }
/****************************************************************/
int Node::getIndex() { return vectorIndex; }
/****************************************************************/
//to be used by printTree in the Tree class to print the leaf info			
void Node::printNode(ostream& out) {
	try{
		out << name;
		
		//there is a branch length
		if (branchLength != -1) { 
			out << ":" << setprecision(4) << branchLength; 
		}
		
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
