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
	KmerAlign(int, int);
	~KmerAlign();
    void align(string, string);
	void align(string, string, vector<int>, vector<int>);
	
private:
    int kmerSize;
    int numKmers;
    int maxKmer;
    Kmer kmerLibrary;
    double calcProb(string A, string B, int overlap, vector<int>, vector<int>);
    double qual_score[47] = {
        -2, -1.58147, -0.996843, -0.695524, -0.507676, -0.38013, -0.289268, -0.222552, -0.172557, -0.134552, -0.105361, -0.0827653, -0.0651742, -0.0514183, -0.0406248, -0.0321336, -0.0254397, -0.0201544, -0.0159759, -0.0126692, -0.0100503, -0.007975, -0.00632956, -0.00502447, -0.00398902, -0.00316729, -0.00251505, -0.00199726, -0.00158615, -0.00125972, -0.0010005, -0.000794644, -0.000631156, -0.000501313, -0.000398186, -0.000316278, -0.00025122, -0.000199546, -0.000158502, -0.0001259, -0.000100005, -7.9436e-05, -6.30977e-05, -5.012e-05, -3.98115e-05, -3.16233e-05, -2.51192e-05};
};

/**************************************************************************************************/

#endif
