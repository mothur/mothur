//
//  optirefmatrix.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef optirefmatrix_hpp
#define optirefmatrix_hpp

#include "optidata.hpp"
#include "optimatrix.h"
#include "subsample.h"

/* Looking to easily access ref, fit and combined information to compare OTU assignments for the references, the sequences to fit, and the merged reference fit OTUs */

class OptiRefMatrix : public OptiData {
    
public:
    
#ifdef UNIT_TEST
    OptiRefMatrix() : OptiData(0.03) { };
#endif
    
    OptiRefMatrix(string, string, string, string, double, float, string); //distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs, refWeightMethod (options: abundance, none, connectivity)
    OptiRefMatrix(string, string, string, string, double, set<string>); //distfile, distFormat, dupsFile, dupsFormat, cutoff, accnosfile refNames
    OptiRefMatrix(string, string, string, string, double, string, string, string, string, string, string); //refdistfile, refname or refcount, refformat, refdistformat, cutoff, fitdistfile, fitname or fitcount, fitformat, fitdistformat, betweendistfile, betweendistformat - files for reference
    ~OptiRefMatrix(){ }
    
    vector<long long> getTranslatedBins(vector<vector<string> >&, vector< vector<long long> >&);
    OptiData* extractMatrixSubset(set<long long>&);
    OptiData* extractMatrixSubset(set<string>&);
    OptiData* extractRefMatrix();
    void randomizeRefs();
    vector<string> getRefSingletonNames();
    
    long long getNumFitTrueSingletons(); //reads that are true singletons (no valid dists in matrix) and are flagged as fit
    long long getNumFitSingletons() { return numFitSingletons; } //user singletons
    long long getNumDists()    { return (numFitDists+numRefDists+numBetweenDists); } //all distances under cutoff
    long long getNumFitDists() { return numFitDists; } //user distances under cutoff
    long long getNumRefDists() { return numRefDists; } //ref distances under cutoff
    
    ListVector* getFitListSingle();
    
    vector<long long> getRefSeqs(); //every ref seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    vector<long long> getFitSeqs(); //every fit seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    
    long long getNumFitSeqs() { return numFitSeqs; } //only Fit seqs that are in fitdistfile and not singletons
    long long getNumFitClose(long long);
    long long getNumRefClose(long long);
    set<long long> getCloseFitSeqs(long long);
    set<long long> getCloseRefSeqs(long long);
    bool isCloseFit(long long, long long, bool&);

protected:
    
    SubSample subsample; 
    map<long long, long long> weights; //seqeunce index in matrix to weight in chosing as reference
    string method, refWeightMethod;
    bool square;
    //a refSingleton or Fitsingleton may not be a true singleton (no valid dists in matrix), but may be a refSeq with no distances to other refs but distances to fitseqs. a fitsingleton may have dists to refs but no dists to other fitseqs.
    long long numFitDists, numRefDists, numRefSingletons, numFitSingletons, numBetweenDists, numSingletons, numFitSeqs;
    float fitPercent;
    
    int readPhylip(string distFile, bool hasName, map<string, string>& names, map<string, long long>& nameAssignment, map<long long, long long>& singletonIndexSwap);
    int readColumn(string distFile, bool hasName, map<string, string>& names, map<string, long long>& nameAssignment, map<long long, long long>& singletonIndexSwap);
    int readFiles(string, string, string, string, string, string, string, string, string, string, string, string);
    int readFiles(string, string, string, string, set<string>&);
    
    map<long long, long long> readColumnSingletons(vector<bool>& singleton, string distFile, map<string, long long>&);
    map<long long, long long> readPhylipSingletons(vector<bool>& singleton, string distFile, long long&,  map<string, long long>& nameAssignment);
    
    vector<bool> isRef; //same size as closeness, this tells us whether a seq with distances in the matrix is a reference
    vector<bool> isSingleRef; ////same size as singletons, this tells us whether a seq WITHOUT distances in the matrix (singleton) is a reference
    
    void calcCounts();
    

    
};


#endif /* optirefmatrix_hpp */
