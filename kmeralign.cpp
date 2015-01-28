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

void KmerAlign::align(string A, string B){
	try {
        
        int nKmersA = A.length() - kmerSize + 1;
        vector<int> kmerA(maxKmer, 0);
        
        int kmer;
        
        for(int i=0;i<nKmersA;i++){
            kmer = kmerLibrary.getKmerNumber(A, i);
            if(kmer != 0){
                kmerA[kmer] = i;    //keep track of where each kmer was found
            }
            else{
                kmerA[kmer] = -1;   //if it's been seen more than once assign it the -1 index of doom
            }
        }
        
        int nKmersB = B.length() - kmerSize + 1;
        int position = -1;
        
        vector<int> kmerB(nKmersB, 0);
        for(int i=0;i<nKmersB;i++){
            kmer = kmerLibrary.getKmerNumber(B, i);

            if(kmerA[kmer] > 0){
                position = i;
                break;
            }
        }
        
        if(position != nKmersB){
            if(kmerA[kmer]-position > 0){ //add gaps to the start of B
                int blength = B.length();
                int alength = A.length();
                int numGaps = (kmerA[kmer]-position);
                B = string(numGaps, '-') + B;
                for (int i = 0; i < blength; i++) {  BBaseMap[i+numGaps] = i;   }
                for (int i = 0; i < alength; i++) {  ABaseMap[i] = i;           }
            }
            else if(kmerA[kmer]-position < 0){ //add gaps to the start of A
                int blength = B.length();
                int alength = A.length();
                int numGaps = (position-kmerA[kmer]);
                A = string(numGaps, '-') + A;
                for (int i = 0; i < alength; i++) {  ABaseMap[i+numGaps] = i;   }
                for (int i = 0; i < blength; i++) {  BBaseMap[i] = i;           }
            }
        }
        int diff = B.length() - A.length();
        if(diff > 0){
            A = A + string(diff, '-');
        }
        else if(diff < 0){
            B = B + string(-diff, '-');
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


