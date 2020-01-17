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

int TaxonomyNode::getChildIndex(string c){
	map<string, int>::iterator it = children.find(c);
	if(it != children.end())	{	return it->second;			}
	else						{	return -1;					}	
}
/**************************************************************************************************/
