/*
 *  taxonomynode.cpp
 *  
 *
 *  Created by Pat Schloss on 7/8/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

/**************************************************************************************************/

#include "taxonomynode.h"

/**************************************************************************************************/

TaxonomyNode::TaxonomyNode(string n, int l): name(n), level(l){
    m = MothurOut::getInstance();
	parent = -1;
	numChildren = 0;
	numSeqs = 0;
}

/**************************************************************************************************/

void TaxonomyNode::setName(string n)			{	name = n;					}

/**************************************************************************************************/

string TaxonomyNode::getName()					{	return name;				}

/**************************************************************************************************/

void TaxonomyNode::setParent(int p)				{	parent = p;					}

/**************************************************************************************************/

int TaxonomyNode::getParent()					{	return parent;				}

/**************************************************************************************************/

void TaxonomyNode::makeChild(string c, int i)	{	children[c] = i;			}


/**************************************************************************************************/

map<string, int> TaxonomyNode::getChildren()	{	return children;			}

/**************************************************************************************************/

int TaxonomyNode::getChildIndex(string c){
	map<string, int>::iterator it = children.find(c);
	if(it != children.end())	{	return it->second;			}
	else						{	return -1;					}	
}

/**************************************************************************************************/

int	TaxonomyNode::getNumKids()					{	return (int)children.size();		}

/**************************************************************************************************/

int	TaxonomyNode::getNumSeqs()					{	return numSeqs;				}

/**************************************************************************************************/

void TaxonomyNode::setTotalSeqs(int n)			{	totalSeqs = n;				}

/**************************************************************************************************/

int TaxonomyNode::getLevel()					{	return level;				}

/**************************************************************************************************/
