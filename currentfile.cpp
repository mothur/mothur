/*
 *  currentfile.cpp
 *  Mothur
 *
 *  Created by westcott on 3/15/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "currentfile.h"

/******************************************************/
CurrentFile* CurrentFile::getInstance() {
	if(_uniqueInstance == 0) {
		_uniqueInstance = new CurrentFile();
	}
	return _uniqueInstance;
}
/*********************************************************************************************/
CurrentFile::CurrentFile() { 
	m = MothurOut::getInstance();
	
	phylipfile = "";
	columnfile = "";
	listfile = "";
	rabundfile = "";
	sabundfile = "";
	namefile = "";
	groupfile = "";
	designfile = "";
	orderfile = "";
	treefile = "";
	sharedfile = "";
	ordergroupfile = "";
	relabundfile = "";
	fastafile = "";
	qualfile = "";
	sfffile = "";
	oligosfile = "";
}
/*********************************************************************************************/
CurrentFile::~CurrentFile() { 
	_uniqueInstance = 0; 
}
/*********************************************************************************************/


