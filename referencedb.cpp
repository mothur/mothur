/*
 *  referencedb.cpp
 *  Mothur
 *
 *  Created by westcott on 6/29/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "referencedb.h"

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
	setSavedTaxonomy("");
}
/*******************************************************
ReferenceDB::~ReferenceDB() { myInstance = NULL; }
/*******************************************************/

