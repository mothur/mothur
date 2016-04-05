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


/* This class is a parent to bayesian, knn.  */

#include "mothur.h"
#include "database.hpp"
#include "phylotree.h"

class Sequence;

/**************************************************************************************************/

class Classify {

public:
	Classify();
	virtual ~Classify(){};
	virtual string getTaxonomy(Sequence*) = 0;
	virtual string getSimpleTax()  { return simpleTax;	}
	virtual bool getFlipped()  { return flipped;	}
	virtual void generateDatabaseAndNames(string, string, string, int, float, float, float, float);
	virtual void setDistName(string s) {} //for knn, so if distance method is selected with knn you can create the smallest distance file in the right place.
    int getMaxLevel() { return maxLevel; }
	
protected:

	map<string, string> taxonomy;  //name maps to taxonomy
	map<string, int>::iterator itTax;
	map<string, string>::iterator it;
	Database* database;
	PhyloTree* phyloTree;
	
	string taxFile, templateFile, simpleTax;
	vector<string> names;
	int threadID, numLevels, numTaxa, maxLevel;
	bool flip, flipped, shortcuts;
	
	int readTaxonomy(string);
	vector<string> parseTax(string);
    double getLogExpSum(vector<double>, int&);
	MothurOut* m;
	
};

/**************************************************************************************************/

#endif

