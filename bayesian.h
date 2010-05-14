#ifndef BAYESIAN_H
#define BAYESIAN_H

/*
 *  bayesian.h
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "classify.h"

/**************************************************************************************************/

class Bayesian : public Classify {
	
public:
	Bayesian(string, string, string, int, int, int);
	~Bayesian();
	
	string getTaxonomy(Sequence*);
	
private:
	vector< vector<float> > wordGenusProb;	//vector of maps from genus to probability
												//wordGenusProb[0][392] = probability that a sequence within genus that's index in the tree is 392 would contain kmer 0;
	
	vector<int> genusTotals;
	vector<int> genusNodes;  //indexes in phyloTree where genus' are located
	
	int kmerSize, numKmers, confidenceThreshold, iters;
	
	string bootstrapResults(vector<int>, int, int);
	int getMostProbableTaxonomy(vector<int>);
	void readProbFile(ifstream&, ifstream&, string, string);
	
};

/**************************************************************************************************/

#endif

