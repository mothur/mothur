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
	m = MothurOut::getInstance();
	//initialize node
	name = "";
	branchLength = -1;
	parent = -1;
	lchild = -1;
	rchild = -1;
	length2leaf = 0.0;
	label = "";
	
}
/****************************************************************/
void Node::setName(string Name) {  name = Name; }
/****************************************************************/
void Node::setGroup(vector<string> groups)  { group =groups; }
/****************************************************************/
void Node::setBranchLength(float l) { branchLength = l; }
/****************************************************************/
void Node::setLabel(string l) { label = l; }
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
vector<string> Node::getGroup() { return group; }
/****************************************************************/
float Node::getBranchLength() { return branchLength; }
/****************************************************************/
string Node::getLabel() { return label; }
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
		m->mothurOut(name + " " + toString(parent) + " " + toString(lchild) + " " + toString(rchild) + " \n");
		
    }
	catch(exception& e) {
		m->errorOut(e, "Node", "printNode");
		exit(1);
	}
}
/****************************************************************/	
