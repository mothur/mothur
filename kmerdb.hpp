#ifndef KMERDB_HPP
#define KMERDB_HPP

/*
 *  kmerdb.h
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"

class KmerDB : public Database {
	
public:
	KmerDB(string, int);
	Sequence* findClosestSequence(Sequence*);

private:
	void generateKmerDB(string);
	void readKmerDB(string, ifstream&);
	int kmerSize;
	int maxKmer;
	vector<vector<int> > kmerLocations;
};

#endif
