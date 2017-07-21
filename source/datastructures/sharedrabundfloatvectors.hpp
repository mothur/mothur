//
//  sharedrabundfloatvectors.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef sharedrabundfloatvectors_hpp
#define sharedrabundfloatvectors_hpp

#include "datavector.hpp"
#include "rabundvector.hpp"
#include "rabundfloatvector.hpp"
#include "sharedordervector.h"

/*  DataStructure for a relabund file. */
//********************************************************************************************************************
inline bool compareRAbundFloats(RAbundFloatVector* left, RAbundFloatVector* right){ return (left->getGroup() < right->getGroup()); }
//********************************************************************************************************************


class SharedRAbundFloatVectors : public DataVector {
    
public:
    SharedRAbundFloatVectors() : DataVector() { label = ""; numBins = 0; }
    SharedRAbundFloatVectors(ifstream&);
    SharedRAbundFloatVectors(SharedRAbundFloatVectors& bv) : DataVector(bv), numBins(bv.numBins) {
        vector<RAbundFloatVector*> data = bv.getSharedRAbundFloatVectors();
        for (int i = 0; i < data.size(); i++) { push_back(data[i]); }
        eliminateZeroOTUS();
        setLabel(bv.getLabel());
    }
    ~SharedRAbundFloatVectors() { for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear(); }
    
    float getOTUTotal(int bin);
    vector<float> getOTU(int bin);
    void setLabel(string l);
    float removeOTUs(vector<int> bins);
    float removeOTU(int bin);
    float get(int bin, string group);
    void set(int bin, float binSize, string group);
    
    int push_back(RAbundFloatVector*);
    void print(ostream&);
    void printHeaders(ostream&);
    void removeGroups(vector<string> g);
    int removeGroups(int minSize, bool silent=false);  // removes any groups with numSeqs < minSize
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedRAbundFloatVectors.\n"); }
    void clear() { for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear(); groupNames.clear(); numBins = 0; }
    int size() { return lookup.size();  }
    int getNumGroups() { return lookup.size(); }
    int getNumBins() { return numBins; }
    float getNumSeqs(string); //group
    float getNumSeqsSmallestGroup();
    
    vector<RAbundVector*> getSharedRAbundVectors();
    vector<RAbundFloatVector*> getSharedRAbundFloatVectors();
    RAbundVector getRAbundVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*);
    
    vector<string> getNamesGroups();
    vector<int> eliminateZeroOTUS(); //run after push_backs if groups are chosen
    
private:
    vector<RAbundFloatVector*> lookup;
    map<string, int> groupNames;
    
    int numBins;

};

#endif /* sharedrabundfloatvectors_hpp */


