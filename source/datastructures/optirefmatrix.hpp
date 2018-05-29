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

/* Looking to easily access ref, fit and combined information to compare OTU assignments for the references, the sequences to fit, and the merged reference fit OTUs */

class OptiRefMatrix : public OptiData {
    
public:
    
    OptiRefMatrix(string, string, string, string, double, string, string, string, string, string, string); //refdistfile, refname or refcount, refformat, refdistformat, cutoff, fitdistfile, fitname or fitcount, fitformat, fitdistformat, betweendistfile, betweendistformat
    ~OptiRefMatrix(){ }
    
    vector<int> getTranslatedBins(vector<vector<string> >&, vector< vector<int> >&);
    
    long long getNumFitSingletons() { return numFitSingletons; } //user singletons
    long long getNumRefSingletons() { return numRefSingletons; } //reference singletons
    
    long long getNumDists()    { return (numFitDists+numRefDists+numBetweenDists); } //all distances under cutoff
    long long getNumFitDists() { return numFitDists; } //user distances under cutoff
    long long getNumRefDists() { return numRefDists; } //ref distances under cutoff
    
    //ListVector* getRefListSingle();
    ListVector* getFitListSingle();
    
    vector<int> getRefSeqs(); //every ref seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    vector<int> getFitSeqs(); //every fit seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    
    long long getNumUniqueFitSeqs() { return numFitSeqs; } //only Fit seqs that are in fitdistfile
    int getNumFitClose(int);
    int getNumRefClose(int);
    set<int> getCloseFitSeqs(int);
    set<int> getCloseRefSeqs(int);
    bool isCloseFit(int, int, bool&);

protected:
    
    string refdistfile, refnamefile, refcountfile, refformat, refdistformat, fitdistfile, fitnamefile, fitcountfile, fitformat, fitdistformat, betweendistfile, betweendistformat, method;
    bool square;
    long long numFitDists, numRefDists, numRefSingletons, numFitSingletons, numBetweenDists, numSingletons;
    long long length, refEnd, refSingletonsEnd, numFitSeqs;  //refEnd is index of first fit seq
    
    int readPhylip(string distFile, bool hasName, map<string, string>& names, map<string, int>& nameAssignment, map<int, int>& singletonIndexSwap);
    int readColumn(string distFile, bool hasName, map<string, string>& names, map<string, int>& nameAssignment, map<int, int>& singletonIndexSwap);
    int readFiles();
    
    map<int, int> readColumnSingletons(vector<bool>& singleton, string namefile, string countfile, string distFile, int&, map<string, int>&);
    void readColumnSingletons(vector<bool>& singleton, string distFile, map<string, int>& nameAssignment);
    map<int, int> readPhylipSingletons(vector<bool>& singleton, string distFile, int&,  map<string, int>& nameAssignment);
    
};


#endif /* optirefmatrix_hpp */
