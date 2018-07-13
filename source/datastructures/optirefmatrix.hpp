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

/* Looking to easily access ref, fit and combined information to compare OTU assignments for the references, the sequences to fit, and the merged reference fit OTUs */

class OptiRefMatrix : public OptiData {
    
public:
    
    OptiRefMatrix(string, string, string, string, double, float); //distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
    OptiRefMatrix(string, string, string, string, double, string, string, string, string, string, string); //refdistfile, refname or refcount, refformat, refdistformat, cutoff, fitdistfile, fitname or fitcount, fitformat, fitdistformat, betweendistfile, betweendistformat - files for reference
    ~OptiRefMatrix(){ }
    
    vector<int> getTranslatedBins(vector<vector<string> >&, vector< vector<int> >&);
    OptiData* extractMatrixSubset(set<int>&);
    OptiData* extractRefMatrix();
    void randomizeRefs();
    vector<string> getRefSingletonNames();
    
    long long getNumFitTrueSingletons(); //reads that are true singletons (no valid dists in matrix) and are flagged as fit
    long long getNumFitSingletons() { return numFitSingletons; } //user singletons
    long long getNumDists()    { return (numFitDists+numRefDists+numBetweenDists); } //all distances under cutoff
    long long getNumFitDists() { return numFitDists; } //user distances under cutoff
    long long getNumRefDists() { return numRefDists; } //ref distances under cutoff
    
    ListVector* getFitListSingle();
    
    vector<int> getRefSeqs(); //every ref seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    vector<int> getFitSeqs(); //every fit seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    
    long long getNumFitSeqs() { return numFitSeqs; } //only Fit seqs that are in fitdistfile and not singletons
    int getNumFitClose(int);
    int getNumRefClose(int);
    set<int> getCloseFitSeqs(int);
    set<int> getCloseRefSeqs(int);
    bool isCloseFit(int, int, bool&);

protected:
    
    string method;
    bool square;
    //a refSingleton or Fitsingleton may not be a true singleton (no valid dists in matrix), but may be a refSeq with no distances to other refs but distances to fitseqs. a fitsingleton may have dists to refs but no dists to other fitseqs.
    long long numFitDists, numRefDists, numRefSingletons, numFitSingletons, numBetweenDists, numSingletons, numFitSeqs;
    float fitPercent;
    
    int readPhylip(string distFile, bool hasName, map<string, string>& names, map<string, int>& nameAssignment, map<int, int>& singletonIndexSwap);
    int readColumn(string distFile, bool hasName, map<string, string>& names, map<string, int>& nameAssignment, map<int, int>& singletonIndexSwap);
    int readFiles(string, string, string, string, string, string, string, string, string, string, string, string);
    int readFiles(string, string, string, string);
    
    map<int, int> readColumnSingletons(vector<bool>& singleton, string namefile, string countfile, string distFile, int&, map<string, int>&);
    void readColumnSingletons(vector<bool>& singleton, string distFile, map<string, int>& nameAssignment);
    map<int, int> readPhylipSingletons(vector<bool>& singleton, string distFile, int&,  map<string, int>& nameAssignment);
    //void assignReferences(vector<string>);
    
    vector<bool> isRef;
    vector<bool> isSingleRef;
    
};


#endif /* optirefmatrix_hpp */
