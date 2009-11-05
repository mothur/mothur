#ifndef CLASSIFY_H
#define CLASSIFY_H

/*
 *  classify.h
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


/* This class is a parent to phylotyp, bayesian, knn.  */

#include "mothur.h"
#include "database.hpp"



class Sequence;

/**************************************************************************************************/

class Classify {

public:
	Classify(string, string, string, int, int, int, int, int);
	Classify(){};
	
	virtual ~Classify(){};
	virtual string getTaxonomy(Sequence*) = 0;
	
protected:

	map<string, string> taxonomy;  //name maps to taxonomy
	map<string, string>::iterator it;
	Database* database;
	
	string taxFile, templateFile;
	vector<string> names;
	
	void readTaxonomy(string);
		
};

/**************************************************************************************************/

#endif

