#ifndef KMERALIGN_N
#define KMERALIGN_N

/*
 *  kmeralign.h
 *
 *
 *  Created by Pat Schloss on 4/6/14.
 *  Copyright 2014 Patrick D. Schloss. All rights reserved.
 *
 *	This class is an Alignment child class that implements a kmer-based pairwise alignment algorithm
 *  for making contigs of reads without insertions
 *
 *
 */

#include "alignment.hpp"
#include "kmer.hpp"

#        define PHREDMAX 46
#        define PHREDCLAMP(x) ((x) > PHREDMAX ? PHREDMAX : ((x) < 0 ? 0 : (x)))


/**************************************************************************************************/

class KmerAlign : public Alignment {
	
public:
	KmerAlign(int);
	~KmerAlign();
    void align(string, string, bool createBaseMap=false);
	
private:
    int kmerSize;
    int maxKmer;
    Kmer kmerLibrary;
    double calcProb(string A, string B, int overlap);
};

/**************************************************************************************************/

#endif
