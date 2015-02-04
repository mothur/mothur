//
//  kmeralign.cpp
//  Mothur
//
//  Created by Pat Schloss on 4/6/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "kmeralign.h"
#include "kmer.hpp"
#include "alignment.hpp"

/**************************************************************************************************/

KmerAlign::KmerAlign(int k, int nk) : numKmers(nk), kmerSize(k), kmerLibrary(k), Alignment() {
	try {
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        //maxKmer = kmerLibrary.getMaxKmer();
        maxKmer = power4s[kmerSize]+1;// (int)pow(4.,k)+1;
	}
	catch(exception& e) {
		m->errorOut(e, "KmerAlign", "KmerAlign");
		exit(1);
	}
}
/**************************************************************************************************/

KmerAlign::~KmerAlign(){	/*	do nothing	*/	}

/**************************************************************************************************/
//modelled after pandaseqs kmer align, assemble.c
void KmerAlign::align(string A, string B, vector<int> AQual, vector<int> BQual){
	try {
        int aLength = A.length();
        int bLength = B.length();
        int maxOverlap = aLength;
        if (bLength < aLength) { maxOverlap = bLength; }
        maxOverlap -= 2;
        
        int nKmersA = A.length() - kmerSize + 1;
        vector< vector<int> > kmerseen;
        //set all kmers to unseen
        kmerseen.resize(maxKmer);
        for (int i = 0; i < maxKmer; i++) { kmerseen[i].resize(numKmers, 0); }
        
        int kmer;
        /* Scan forward sequence building k-mers and appending the position to kmerseen[k] */
        for(int i=1;i<nKmersA;i++){
            kmer = kmerLibrary.getKmerNumber(A, i);
            cout << i << '\t' << kmer << endl;
            int j;
            for (j = 0; j < numKmers && kmerseen[kmer][j] != 0; j++) ;
            printf("kmerseen: %i, %i, %i\n", numKmers, j, (kmer * numKmers) + j);
            if (j == numKmers) {
                /* If we run out of storage, we lose k-mers. */
                m->mothurOut("[WARNING]: Inssufficent number of kmer instances. Try increasing the kmer parameter.\n");
            } else {
                kmerseen[kmer][j] = i;
            }
        }
        int nKmersB = B.length() - kmerSize + 1;
        
        /* Scan reverse sequence building k-mers. For each position in the forward sequence for this kmer (i.e., kmerseen[k]), flag that we should check the corresponding overlap. */
        set<int> overlaps;
        for(int i=1;i<nKmersB;i++){
            kmer = kmerLibrary.getKmerNumber(B, i);
            printf("kmerseen reverse: %i\n", kmer);
            for (int j = 0; j < numKmers && kmerseen[kmer][j] != 0; j++) { //for as many instances as we saw of this kmer, recoding overlap
                int index = aLength + bLength - (bLength-i) - kmerseen[kmer][j] - 2;
                if (index <= maxOverlap) { overlaps.insert(index); }
                printf("kmerseen reverse: %i, %i, %i, %i\n", index, kmerseen[kmer][j], (bLength-i), i);
            }
        }
        
        //find best overlap
        double bestProb = -1.0;
        int bestOverlap = -1;
        for (set<int>::iterator it = overlaps.begin(); it != overlaps.end(); it++) {
            int overlap = (*it) + 2; //2 = minoverlap
            double probability = calcProb(A, B, overlap, AQual, BQual);
            printf("overlap prob: %i, %i, %f\n", *it, overlap, probability);
            if (probability > bestProb && overlap >= 2) {
                bestProb = probability;
                bestOverlap = overlap;
            }
        }
        cout << bestOverlap << '\t' << (aLength-(bestOverlap)) << endl;
        if(bestOverlap != -1){
            if((aLength-bestOverlap) > 0){ //add gaps to the start of B
                int numGaps = (aLength-bestOverlap);
                B = string(numGaps, '-') + B;
                for (int i = 0; i < bLength; i++) {  BBaseMap[i+numGaps] = i;   }
                for (int i = 0; i < aLength; i++) {  ABaseMap[i] = i;           }
            }
        }
        int diff = B.length() - A.length();
        if(diff > 0){
            A = A + string(diff, '-');
        }
        
        seqAaln = A;
        seqBaln = B;
        pairwiseLength = seqAaln.length();
        
    }
	catch(exception& e) {
		m->errorOut(e, "KmerAlign", "align");
		exit(1);
	}
    
}
/**************************************************************************************************/
//modelled after pandaseqs kmer align, assemble.c
double KmerAlign::calcProb(string A, string B, int overlap, vector<int> forwardQual, vector<int> reverseQual){
    try {
        double prob = 0;
        int aLength = A.length();
        int bLength = B.length();
        double randomBase = log(0.25);
        
        for (int i = 0; i < overlap; i++) {
            int findex = aLength + i - overlap;
            int rindex = bLength - (bLength-i) - 1;
            if (findex < 0 || rindex < 0 || findex >= aLength || rindex >= bLength)
                continue;
            char f = A[findex];
            char r = B[rindex];
            if ((f == 'N') || (r == 'N')) {
                prob -= randomBase;
            } else if ((f & r) != 0) {
                //probability += qual_match_pear[PHREDCLAMP(forward[findex].qual)][PHREDCLAMP(forward[rindex].qual)];
                double p = qual_score[PHREDCLAMP(forwardQual[findex])];
                double q = qual_score[PHREDCLAMP(forwardQual[rindex])];
                prob += (1.0 - (1.0 - q) * p / 3.0 - (1.0 - p) * q / 3.0 - 2.0 * (1.0 - p) * (1.0 - q) / 9.0);
            } else {
                double p = qual_score[PHREDCLAMP(reverseQual[findex])];
                double q = qual_score[PHREDCLAMP(reverseQual[rindex])];
                prob += (1 - p) * q / 3 + (1 - q) * p / 3 + p * q / 2;
            }
        }
        
        return prob;
    }
    catch(exception& e) {
        m->errorOut(e, "KmerAlign", "calcProb");
        exit(1);
    }
}
/**************************************************************************************************/

void KmerAlign::align(string A, string B){
    try {
        m->mothurOut("[ERROR]: kmerAlign requires quality scores.\n");
        m->control_pressed = true;
    }
    catch(exception& e) {
        m->errorOut(e, "KmerAlign", "align");
        exit(1);
    }
    
}

/**************************************************************************************************/


