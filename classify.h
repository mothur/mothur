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
	Classify(string, string, string, int, float, float, float, float);
	
	virtual ~Classify(){  delete phyloTree; delete database;  };
	virtual string getTaxonomy(Sequence*) = 0;
	//virtual map<string, int> getConfidenceScores() { return taxConfidenceScore; }
	//virtual vector<string> parseTax(string);
	virtual string getSimpleTax()  { return simpleTax;	}
	
protected:

	map<string, string> taxonomy;  //name maps to taxonomy
	//map<string, int> genusCount;  //maps genus to count - in essence a list of how many seqs are in each taxonomy
	map<string, int>::iterator itTax;
	map<string, string>::iterator it;
	Database* database;
	PhyloTree* phyloTree;
	
	string taxFile, templateFile, simpleTax;
	vector<string> names;
	
	void readTaxonomy(string);
	vector<string> parseTax(string);
	MothurOut* m;
};

/**************************************************************************************************/

#endif

