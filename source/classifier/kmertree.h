//
//  kmerTree.h
//  pdsBayesian
//
//  Created by Patrick Schloss on 4/3/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pdsBayesian_kmerTree_h
#define pdsBayesian_kmerTree_h

#include "classify.h"

class KmerNode;

class KmerTree : public Classify {
	
public:
	KmerTree(string, string, int, int);
	~KmerTree();
	
    string getTaxonomy(Sequence*, string&, bool&);

private:
    int addTaxonomyToTree(string, string, vector<int>&);
	vector<int> ripKmerProfile(string);
	int getMinRiskIndexKmer(vector<int>&, vector<int>&, vector<double>&);
	int aggregateThetas();
	int sanityCheck(vector<vector<int> >&, vector<int>&);

	int kmerSize;
	int numPossibleKmers, confidenceThreshold;
	vector<KmerNode*> tree;

};

#endif
