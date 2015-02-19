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

KmerAlign::KmerAlign(int k) : kmerSize(k), kmerLibrary(k), Alignment() {
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
void KmerAlign::align(string A, string B){
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
        //for (int i = 0; i < maxKmer; i++) { kmerseen[i].resize(numKmers, 0); }
        
        int kmer;
        /* Scan forward sequence building k-mers and appending the position to kmerseen[k] */
        for(int i=0;i<nKmersA;i++){
            kmer = kmerLibrary.getKmerNumber(A, i);
            if (kmer != (maxKmer-1)) {  kmerseen[kmer].push_back(i);  }
        }
        int nKmersB = B.length() - kmerSize + 1;
        
        /* Scan reverse sequence building k-mers. For each position in the forward sequence for this kmer (i.e., kmerseen[k]), flag that we should check the corresponding overlap. */
        set<int> overlaps;
        for(int i=0;i<nKmersB;i++){
            kmer = kmerLibrary.getKmerNumber(B, i);
            for (int j = 0; j < kmerseen[kmer].size(); j++) {  //for as many instances as we saw of this kmer, recoding overlap
                int index = aLength + bLength - (bLength-i) - kmerseen[kmer][j] - 2;
                if (index <= maxOverlap) { overlaps.insert(index); }
            }
        }
        
        if (overlaps.size() == 0) { //add all possible overlaps
            for (int i = 0; i <= maxOverlap; i++) {  overlaps.insert(i);  } //minoverlap=2
        }
        
        //find best overlap
        double bestProb = -1.38629 * (aLength+bLength);
        int bestOverlap = -1;
        for (set<int>::iterator it = overlaps.begin(); it != overlaps.end(); it++) {
            int index = *it;
            int overlap = index + 2; //2 = minoverlap
            double probability = calcProb(A, B, overlap);
            //printf("overlap prob: %i, %f\n", overlap, probability);
            if (probability > bestProb && overlap >= 2) {
                bestProb = probability;
                bestOverlap = overlap;
            }
        }
        //printf("best overlap prob: %i, %f\n", bestOverlap, bestProb);
        if(bestOverlap != -1){
            if((aLength-bestOverlap) > 0){ //add gaps to the start of B
                int numGaps = (aLength-bestOverlap);
                B = string(numGaps, '-') + B;
                for (int i = 0; i < bLength; i++) {  BBaseMap[i+numGaps] = i;   }
                for (int i = 0; i < aLength; i++) {  ABaseMap[i] = i;           }
            }else {
                for (int i = 0; i < bLength; i++) {  BBaseMap[i] = i;   }
                for (int i = 0; i < aLength; i++) {  ABaseMap[i] = i;   }
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
double KmerAlign::calcProb(string A, string B, int overlap){
    try {
        double prob = 0;
        int aLength = A.length();
        int bLength = B.length();
        int unknown, match, mismatch; unknown = 0; match = 0; mismatch = 0;
    
        for (int i = 0; i < overlap; i++) {
            int findex = aLength + i - overlap;
            int rindex = i;
            if (findex < 0 || rindex < 0 || findex >= aLength || rindex >= bLength)
                continue;
            char f = A[findex];
            char r = B[rindex];
            if ((f == 'N') || (r == 'N')) {
                unknown++;
            } else if (r == f) {
                match++;
            } else {
                mismatch++;
            }
        }
        //ln(0.25 * (1 - 2 * 0.36 + 0.36 * 0.36))
        double pmatch = -2.278869;
        //ln((3 * 0.36 - 2 * 0.36 * 0.36) / 18.0)
        double pmismatch = -3.087848;
        if (overlap >= aLength && overlap >= bLength) {
            prob = (-1.38629 * unknown + match * pmatch + mismatch * pmismatch);
        } else {
            prob = (-1.38629 * (aLength + bLength - 2 * overlap + unknown) + match * pmatch + mismatch * pmismatch);
        }

        
        return prob;
    }
    catch(exception& e) {
        m->errorOut(e, "KmerAlign", "calcProb");
        exit(1);
    }
}
/**************************************************************************************************/


