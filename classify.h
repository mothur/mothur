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
	
	virtual ~Classify(){  delete phyloTree; if (database != NULL) {  delete database; } };
	virtual string getTaxonomy(Sequence*) = 0;
	//virtual map<string, int> getConfidenceScores() { return taxConfidenceScore; }
	//virtual vector<string> parseTax(string);
	virtual string getSimpleTax()  { return simpleTax;	}
	virtual void generateDatabaseAndNames(string, string, string, int, float, float, float, float);
	
protected:

	map<string, string> taxonomy;  //name maps to taxonomy
	//map<string, int> genusCount;  //maps genus to count - in essence a list of how many seqs are in each taxonomy
	map<string, int>::iterator itTax;
	map<string, string>::iterator it;
	Database* database;
	PhyloTree* phyloTree;
	
	string taxFile, templateFile, simpleTax;
	vector<string> names;
	
	int readTaxonomy(string);
	vector<string> parseTax(string);
	MothurOut* m;
	
};

/**************************************************************************************************/

#endif

