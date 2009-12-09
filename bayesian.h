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
	Bayesian(string, string, string, int, int, bool);
	~Bayesian() {};
	
	string getTaxonomy(Sequence*);
	//map<string, int> getConfidenceScores() { return taxConfidenceScore; }
	
private:
	//map<string, vector<double> > taxonomyProbability;	//probability of a word being in a given taxonomy. 
													//taxonomyProbability[bacteria;][0] = probabtility that a sequence within bacteria; would contain kmer 0;
													//taxonomyProbability[bacteria;][1] = probabtility that a sequence within bacteria; would contain kmer 1...
	
	vector< vector<float> > wordGenusProb;	//vector of maps from genus to probability
												//wordGenusProb[0][392] = probability that a sequence within genus that's index in the tree is 392 would contain kmer 0;
	
	vector<int> genusTotals;
	vector<int> genusNodes;  //indexes in phyloTree where genus' are located
	
	int kmerSize, numKmers, confidenceThreshold;
	bool probs;
	
	string bootstrapResults(vector<int>, int, int);
	int getMostProbableTaxonomy(vector<int>);
	void readProbFile(ifstream&, ifstream&);
	//map<string, int> parseTaxMap(string);
	
};

/**************************************************************************************************/

#endif

