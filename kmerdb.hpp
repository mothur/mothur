#ifndef KMERDB_HPP
#define KMERDB_HPP

/*
 *  kmerdb.h
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class is a child class of the Database class, which stores the template sequences as a kmer table and provides
 *	a method of searching the kmer table for the sequence with the most kmers in common with a query sequence.
 *	kmerLocations is the primary storage variable that is a two-dimensional vector where each row represents the
 *	different number of kmers and each column contains the index to sequences that use that kmer.
 *
 *	Construction of an object of this type will first look for an appropriately named database file and if it is found
 *	then will read in the database file (readKmerDB), otherwise it will generate one and store the data in memory
 *	(generateKmerDB)

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
