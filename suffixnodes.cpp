/*
 *  SuffixNodes.cpp
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

#include "suffixnodes.hpp"


//********************************************************************************************************************

inline char deCodeSequence(char code){
	
	if(code == '0')			{	return 'a';	}	//	this method allows us to go from the int string to a char string;
	else if(code == '1')	{	return 'c';	}	//	it's only really useful if we want to print out the tree
	else if(code == '2')	{	return 'g';	}
	else if(code == '3')	{	return 't';	}
	else if(code == '4')	{	return 'n';	}
	else					{	return '$';	}
	
}

//********************************************************************************************************************

SuffixNode::SuffixNode(int parent, int start, int end) : 
		parentNode(parent),			//	we store the parent node as an int
		startCharPosition(start),	//	the suffix tree class will hold the sequence that the startCharPosition and 
		endCharPosition(end)		//	endCharPosition indices correspond to
		{	/*	do nothing	*/			}


void SuffixNode::setChildren(char, int)			{	/*	do nothing	*/			}	//	there's no children in a leaf
int SuffixNode::getNumChildren()				{	return 0;					}	//	ditto
void SuffixNode::eraseChild(char)				{	/*	do nothing	*/			}	//	ditto
int SuffixNode::getChild(char)					{	return -1;					}	//	ditto
void SuffixNode::setSuffixNode(int)				{	/*	do nothing	*/			}	//	there's no suffix node in a leaf
int SuffixNode::getSuffixNode()					{	return -1;					}	//	ditto
int SuffixNode::getParentNode()					{	return parentNode;			}
void SuffixNode::setParentNode(int number)		{	parentNode = number;		}
int SuffixNode::getStartCharPos()				{	return startCharPosition;	}
void SuffixNode::setStartCharPos(int start)		{	startCharPosition = start;	}
int SuffixNode::getEndCharPos()					{	return endCharPosition;		}	

//********************************************************************************************************************

SuffixLeaf::SuffixLeaf(int parent, int start, int end) : SuffixNode(parent, start, end) {	/*	do nothing	*/	}


void SuffixLeaf::print(string sequence, int nodeNumber){
	
	cout << this << '\t' << parentNode << '\t' << nodeNumber << '\t' <<
	-1 << '\t' << startCharPosition << '\t' << endCharPosition << '\t';
	
	cout << '\'';
	for(int i=startCharPosition;i<=endCharPosition;i++){
		cout << deCodeSequence(sequence[i]);
	}
	cout << '\'' << endl;
}

//********************************************************************************************************************

SuffixBranch::SuffixBranch(int parent, int start, int end) : SuffixNode(parent, start, end), suffixNode(-1){
		childNodes.assign(6, -1);
}
	
void SuffixBranch::print(string sequence, int nodeNumber){						//	this method is different that than
	cout << this << '\t' << parentNode << '\t' << nodeNumber << '\t' <<			//	of a leaf because it prints out a
	suffixNode << '\t' << startCharPosition << '\t' << endCharPosition << '\t';	//	value for the suffix node	
	
	cout << '\'';
	for(int i=startCharPosition;i<=endCharPosition;i++){
		cout << deCodeSequence(sequence[i]);
	}
	cout << '\'' << endl;
}

//	we can access the children by subtracting '0' from the the char value from the string, the difference is an int
//	value and the index we need to access.
void SuffixBranch::eraseChild(char base)	{	childNodes[base - '0'] = -1;	}	//to erase set the child index to -1
void SuffixBranch::setChildren(char base, int nodeIndex){	childNodes[base - '0'] = nodeIndex;	}
void SuffixBranch::setSuffixNode(int nodeIndex){	suffixNode = nodeIndex;		}
int SuffixBranch::getSuffixNode()			{	return suffixNode;				}
int SuffixBranch::getChild(char base)		{	return childNodes[base - '0'];	}
	
//********************************************************************************************************************
