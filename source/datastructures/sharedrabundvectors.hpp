//
//  sharedrabundvectors.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef sharedrabundvectors_hpp
#define sharedrabundvectors_hpp

#include "datavector.hpp"
#include "rabundvector.hpp"
#include "rabundfloatvector.hpp"
#include "sharedordervector.h"

/*  DataStructure for a shared file. */

//********************************************************************************************************************
inline bool compareRAbunds(RAbundVector* left, RAbundVector* right){ return (left->getGroup() < right->getGroup()); }
//********************************************************************************************************************

class SharedRAbundVectors : public DataVector {
    
public:
    SharedRAbundVectors() : DataVector() {  label = ""; numBins = 0; }
    SharedRAbundVectors(ifstream&);
    SharedRAbundVectors(SharedRAbundVectors& bv) : DataVector(bv), numBins(bv.numBins) {
        vector<RAbundVector*> data = bv.getSharedRAbundVectors();
        for (int i = 0; i < data.size(); i++) { push_back(data[i]); }
        eliminateZeroOTUS();
        setLabel(bv.getLabel());
    }
    ~SharedRAbundVectors() { clear(); }
    
    void setLabel(string l);
    int getOTUTotal(int bin);
    vector<int> getOTU(int bin);
    int push_back(vector<int>, string binLabel=""); //add otu. mothur assumes abunds are in same order as groups.
    int get(int bin, string group);
    void set(int bin, int binSize, string group);
    int push_back(RAbundVector*);
    
    vector<int> eliminateZeroOTUS(); //run after push_backs if groups are chosen
    int removeOTU(int bin);
    int removeOTUs(vector<int> bins);
    void removeGroups(vector<string> g);
    int removeGroups(int minSize, bool silent=false);  // removes any groups with numSeqs < minSize
    int size() { return lookup.size(); }
    int getNumBins() { return numBins; }
    int getNumSeqsSmallestGroup();
    vector<string> getNamesGroups(); //same order as Rabunds
    int getNumGroups() { return lookup.size(); }
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedRAbundVectors.\n"); m->control_pressed = true; }
    void clear() { for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear(); groupNames.clear(); numBins = 0; }
    
    void print(ostream&);
    void printHeaders(ostream&);
    int getNumSeqs(string); //group
    
    RAbundVector getRAbundVector();
    RAbundVector getRAbundVector(string); //group you want the rabund for
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*);
    SharedOrderVector getSharedOrderVector();
    vector<RAbundVector*> getSharedRAbundVectors();
    vector<RAbundFloatVector*> getSharedRAbundFloatVectors();
    
private:
    vector<RAbundVector*> lookup;
    map<string, int> groupNames;
    int numBins;
    
};


#endif /* sharedrabundvectors_hpp */
