/*
 *  suffixtree.cpp
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

#include "sequence.hpp"
#include "suffixnodes.hpp"
#include "suffixtree.hpp"


//********************************************************************************************************************

inline bool compareParents(SuffixNode* left, SuffixNode* right){//	this is necessary to print the tree and to sort the
	return (left->getParentNode() < right->getParentNode());	//	nodes in order of their parent
}

//********************************************************************************************************************

SuffixTree::SuffixTree(){}

//********************************************************************************************************************

SuffixTree::~SuffixTree(){
	for(int i=0;i<nodeVector.size();i++){	delete nodeVector[i];	}	
	nodeVector.clear();
}

//********************************************************************************************************************

void SuffixTree::loadSequence(Sequence seq){
	nodeCounter = 0;							//	initially there are 0 nodes in the tree
	activeStartPosition = 0;
	activeEndPosition = -1;						
	seqName = seq.getName();
	sequence = seq.convert2ints();
	sequence += '5';							//	this essentially concatenates a '$' to the end of the sequence to
	int seqLength = sequence.length();			//	make it a cononical suffix tree
	
	nodeVector.push_back(new SuffixBranch(-1, 0, -1));	//	enter the root of the suffix tree
	
	activeNode = root = 0;
	string hold;
	for(int i=0;i<seqLength;i++){
		addPrefix(i);							//	step through the sequence adding each prefix
	}
}

//********************************************************************************************************************

string SuffixTree::getSeqName()	{
	return seqName;		
}

//********************************************************************************************************************

void SuffixTree::print(){
	vector<SuffixNode*> hold = nodeVector;
	sort(hold.begin(), hold.end(), compareParents);
	cout << "Address\t\tParent\tNode\tSuffix\tStartC\tEndC\tSuffix" << endl;
	for(int i=1;i<=nodeCounter;i++){
		hold[i]->print(sequence, i);
	}
}

//********************************************************************************************************************

int SuffixTree::countSuffixes(string compareSequence, int& minValue){	//	here we count the number of suffix parts 
															//	we need to rewrite a user supplied sequence.  if the 
	int numSuffixes = 0;									//	count exceeds the supplied minValue, bail out.  The
	int seqLength = compareSequence.length();				//	time complexity should be O(L)
	int position = 0;
	
	int presentNode = 0;
	
	while(position < seqLength){		//	while the position in the query sequence isn't at the end...
		
		if(numSuffixes > minValue)	{	return 1000000;		}	//	bail if the count gets too high
		
		int newNode = nodeVector[presentNode]->getChild(compareSequence[position]);	//	see if the current node has a
																//	child that matches the next character in the query
		if(newNode == -1){										
			if(presentNode == 0){	position++;		}			//	if not, go back to the root and increase the count
			numSuffixes++;										//	by one.
			presentNode = 0;
		}
		else{													//	if there is, move to that node and see how far down
			presentNode = newNode;								//	it we can get
			
			for(int i=nodeVector[newNode]->getStartCharPos(); i<=nodeVector[newNode]->getEndCharPos(); i++){
				if(compareSequence[position] == sequence[i]){
					position++;									//	as long as the query and branch agree, keep going
				}
				else{
					numSuffixes++;								//	if there is a mismatch, increase the number of 
					presentNode = 0;							//	suffixes and go back to the root
					break;
				}
			}
		}
		//	if we get all the way through the node we'll go to the top of the while loop and find the child node
		//	that corresponds to what we are interested in		
	}
	numSuffixes--;												//	the method puts an extra count on numSuffixes
	
	if(numSuffixes < minValue)	{	minValue = numSuffixes;	}	//	if the count is less than the previous minValue,
	return numSuffixes;											//	change the value and return the number of suffixes
	
}

//********************************************************************************************************************

void SuffixTree::canonize(){	//	if you have to ask how this works, you don't really want to know and this really
								//	isn't the place to ask.
	if ( isExplicit() == 0 ) {	//	if the node has no children...
		
		int tempNodeIndex = nodeVector[activeNode]->getChild(sequence[activeStartPosition]);
		SuffixNode* tempNode = nodeVector[tempNodeIndex];
		
		int span = tempNode->getEndCharPos() - tempNode->getStartCharPos();
		
		while ( span <= ( activeEndPosition - activeStartPosition ) ) {
			
            activeStartPosition = activeStartPosition + span + 1;
			
			activeNode = tempNodeIndex;
			
            if ( activeStartPosition <= activeEndPosition ) {
				tempNodeIndex = nodeVector[tempNodeIndex]->getChild(sequence[activeStartPosition]);
				tempNode = nodeVector[tempNodeIndex];
				span = tempNode->getEndCharPos() - tempNode->getStartCharPos();
            }
			
        }
    }
}

//********************************************************************************************************************

int SuffixTree::split(int nodeIndex, int position){	//	leaves stay leaves, etc, to split a leaf we make a new interior 
													//	node and reconnect everything
	SuffixNode* node = nodeVector[nodeIndex];					//	get the node that needs to be split
	SuffixNode* parentNode = nodeVector[node->getParentNode()];	//	get it's parent node
	
	parentNode->eraseChild(sequence[node->getStartCharPos()]);	//	erase the present node from the registry of its parent
	
	nodeCounter++;
	SuffixNode* newNode = new SuffixBranch(node->getParentNode(), node->getStartCharPos(), node->getStartCharPos() + activeEndPosition - activeStartPosition);	//	create a new node that will link the parent with the old child
	parentNode->setChildren(sequence[newNode->getStartCharPos()], nodeCounter);//	give the parent the new child
	nodeVector.push_back(newNode);
	
	node->setParentNode(nodeCounter);	//	give the original node the new node as its parent
	newNode->setChildren(sequence[node->getStartCharPos() + activeEndPosition - activeStartPosition + 1], nodeIndex);
	//	put the original node in the registry of the new node's children
	newNode->setSuffixNode(activeNode);//link the new node with the old active node
	
	//	recalculate the startCharPosition of the outermost node
	node->setStartCharPos(node->getStartCharPos() + activeEndPosition - activeStartPosition + 1 );
	
	return node->getParentNode();
}

//********************************************************************************************************************

void SuffixTree::makeSuffixLink(int& previous, int present){
	
//	here we link the nodes that are suffixes of one another to rapidly speed through the tree
	if ( previous > 0 ) {	nodeVector[previous]->setSuffixNode(present);	}
	else				{	/*	do nothing								*/	}
	
    previous = present;
}

//********************************************************************************************************************

void SuffixTree::addPrefix(int prefixPosition){
	
	int lastParentNode = -1;	//	we need to place a new prefix in the suffix tree
	int parentNode = 0;
	
	while(1){
		
		parentNode = activeNode;
		
		if(isExplicit() == 1){	//	if the node is explicit (has kids), try to follow it down the branch if its there...
			if(nodeVector[activeNode]->getChild(sequence[prefixPosition]) != -1){	//	break out and get next prefix...
				break;												
			}
			else{				//	...otherwise continue, we'll need to make a new node later on...
			}
		}
		else{					//	if it's not explicit (no kids), read through and see if all of the chars agree...
			int tempNode = nodeVector[activeNode]->getChild(sequence[activeStartPosition]);
			int span = activeEndPosition - activeStartPosition;
			
			if(sequence[nodeVector[tempNode]->getStartCharPos() + span + 1] == sequence[prefixPosition] ){
				break;			//	if the existing suffix agrees with the new one, grab a new prefix...
			}
			else{
				parentNode = split(tempNode, prefixPosition);	//	... otherwise we need to split the node
			}
			
		}
		
		nodeCounter++;	//	we need to generate a new node here if the kid didn't exist, or we split a node
		SuffixNode* newSuffixLeaf = new SuffixLeaf(parentNode, prefixPosition, sequence.length()-1);
		nodeVector[parentNode]->setChildren(sequence[prefixPosition], nodeCounter);
		nodeVector.push_back(newSuffixLeaf);
		
		makeSuffixLink( lastParentNode, parentNode );		//	make a suffix link for the parent node
		
		if(nodeVector[activeNode]->getParentNode() == -1){	//	move along the start position for the tree
            activeStartPosition++;
        } 
		else {
            activeNode = nodeVector[activeNode]->getSuffixNode();
		}
		canonize();											//	frankly, i'm not entirely clear on what canonize does.
	}
	
	makeSuffixLink( lastParentNode, parentNode );
	activeEndPosition++;									//	move along the end position for the tree
	
	canonize();												//	frankly, i'm not entirely clear on what canonize does.
	
}

//********************************************************************************************************************

