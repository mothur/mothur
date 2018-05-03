/*
 *  kmerNode.cpp
 *  bayesian
 *
 *  Created by Pat Schloss on 10/11/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "kmernode.h"


/**********************************************************************************************************************/

KmerNode::KmerNode(string s, int l, int n) : TaxonomyNode(s, l), kmerSize(n) {
	try {
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        
        numPossibleKmers = power4s[kmerSize];
        numUniqueKmers = 0;
        
        kmerVector.assign(numPossibleKmers, 0);
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "KmerNode");
		exit(1);
	}
}

/**********************************************************************************************************************/

void KmerNode::loadSequence(vector<int>& kmerProfile){
	try {
        for(int i=0;i<numPossibleKmers;i++){
            if (m->getControl_pressed()) { break; }
            if(kmerVector[i] == 0 && kmerProfile[i] != 0)	{	numUniqueKmers++;	}
            
            kmerVector[i] += kmerProfile[i];
        }
        
        numSeqs++;
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "loadSequence");
		exit(1);
	}
}	

/**********************************************************************************************************************/

string KmerNode::getKmerBases(int kmerNumber){
	try {
        //	Here we convert the kmer number into the kmer in terms of bases.
        //
        //	Example:	Score = 915 (for a 6-mer)
        //				Base6 = (915 / 4^0) % 4 = 915 % 4 = 3 => T	[T]
        //				Base5 = (915 / 4^1) % 4 = 228 % 4 = 0 => A	[AT]
        //				Base4 = (915 / 4^2) % 4 = 57 % 4 = 1 => C	[CAT]
        //				Base3 = (915 / 4^3) % 4 = 14 % 4 = 2 => G	[GCAT]
        //				Base2 = (915 / 4^4) % 4 = 3 % 4 = 3 => T	[TGCAT]
        //				Base1 = (915 / 4^5) % 4 = 0 % 4 = 0 => A	[ATGCAT] -> this checks out with the previous method
        
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        
        string kmer = "";
        
        if(kmerNumber == power4s[kmerSize]){//pow(4.,7)){	//	if the kmer number is the same as the maxKmer then it must
            for(int i=0;i<kmerSize;i++){					//	have had an N in it and so we'll just call it N x kmerSize
                kmer += 'N';
            }
        }
        else{
            for(int i=0;i<kmerSize;i++){
                if (m->getControl_pressed()) { return kmer; }
                int nt = (int)(kmerNumber / (float)power4s[i]) % 4;		//	the '%' operator returns the remainder 
                if(nt == 0)		{	kmer = 'A' + kmer;	}				//	from int-based division ]
                else if(nt == 1){	kmer = 'C' + kmer;	}
                else if(nt == 2){	kmer = 'G' + kmer;	}
                else if(nt == 3){	kmer = 'T' + kmer;	}
            }
        }
        return kmer;
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "getKmerBases");
		exit(1);
	}
}

/**************************************************************************************************/

void KmerNode::addThetas(vector<int> newTheta, int newNumSeqs){
	try {
        for(int i=0;i<numPossibleKmers;i++){
            if (m->getControl_pressed()) { break; }
            kmerVector[i] += newTheta[i];		
        }
        
        //	if(alignLength == 0){
        //		alignLength = (int)newTheta.size();
        //		theta.resize(alignLength);
        //		columnCounts.resize(alignLength);
        //	}
        //	
        //	for(int i=0;i<alignLength;i++){	
        //		theta[i].A += newTheta[i].A;		columnCounts[i] += newTheta[i].A;
        //		theta[i].T += newTheta[i].T;		columnCounts[i] += newTheta[i].T;
        //		theta[i].G += newTheta[i].G;		columnCounts[i] += newTheta[i].G;
        //		theta[i].C += newTheta[i].C;		columnCounts[i] += newTheta[i].C;
        //		theta[i].gap += newTheta[i].gap;	columnCounts[i] += newTheta[i].gap;
        //	}
        
        numSeqs += newNumSeqs;
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "addThetas");
		exit(1);
	}
}

/**********************************************************************************************************************/

int KmerNode::getNumUniqueKmers(){
    try {
        if(numUniqueKmers == 0){
            
            for(int i=0;i<numPossibleKmers;i++){
                if (m->getControl_pressed()) { return numUniqueKmers; }
                if(kmerVector[i] != 0){
                    numUniqueKmers++;
                }
                
            }
            
        }
        
        return numUniqueKmers;	
        
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "getNumUniqueKmers");
		exit(1);
	}
}

/**********************************************************************************************************************/

void KmerNode::printTheta(){
	try {
        m->mothurOut(name + "\n");
        for(int i=0;i<numPossibleKmers;i++){
            if(kmerVector[i] != 0){
                m->mothurOut(getKmerBases(i) + '\t' + toString(kmerVector[i]) + "\n");
            }
        }
        m->mothurOutEndLine();	
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "printTheta");
		exit(1);
	}
	
}
/**************************************************************************************************/

double KmerNode::getSimToConsensus(vector<int>& queryKmerProfile){
	try {
        double present = 0;
        
        for(int i=0;i<numPossibleKmers;i++){
            if (m->getControl_pressed()) { return present; }
            if(queryKmerProfile[i] != 0 && kmerVector[i] != 0){
                present++;
            }
        }	
        
        return (present / double(queryKmerProfile.size() - kmerSize + 1));
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "getSimToConsensus");
		exit(1);
	}
}

/**********************************************************************************************************************/

double KmerNode::getPxGivenkj_D_j(vector<int>& queryKmerProfile)	{	
	try {
        double sumLogProb = 0.0000;
        double alpha = 1.0 / (double)totalSeqs;	//flat prior
        //	double alpha = pow((1.0 / (double)numUniqueKmers), numSeqs)+0.0001;	//non-flat prior
        
        for(int i=0;i<numPossibleKmers;i++){
            if (m->getControl_pressed()) { return sumLogProb; }
            if(queryKmerProfile[i] != 0){		//numUniqueKmers needs to be the value from Root;
                sumLogProb += log((kmerVector[i] + alpha) / (numSeqs + numUniqueKmers * alpha));
            }
            
        }
        return sumLogProb;
    }
	catch(exception& e) {
		m->errorOut(e, "KmerNode", "getPxGivenkj_D_j");
		exit(1);
	}

}

/**********************************************************************************************************************/
