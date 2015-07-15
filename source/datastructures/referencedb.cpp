/*
 *  referencedb.cpp
 *  Mothur
 *
 *  Created by westcott on 6/29/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "referencedb.h"

//needed for testing project
//ReferenceDB* ReferenceDB::myInstance;

/******************************************************/
ReferenceDB* ReferenceDB::getInstance()  {
	 if(myInstance == NULL) {
		myInstance = new ReferenceDB();
	 }
	 return myInstance;
 }
/******************************************************/
void ReferenceDB::clearMemory()  {
	referenceSeqs.clear();	
	setSavedReference("");
	for(int i = 0; i < wordGenusProb.size(); i++) { wordGenusProb[i].clear(); }
	wordGenusProb.clear();
	WordPairDiffArr.clear();
	setSavedTaxonomy("");
}
/*******************************************************
ReferenceDB::~ReferenceDB() { myInstance = NULL; }
*******************************************************/

