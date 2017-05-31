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

class SharedRAbundVectors : public DataVector {
    
public:
    SharedRAbundVectors() : DataVector() {  label = ""; numBins = 0; }
    SharedRAbundVectors(ifstream&);
    SharedRAbundVectors(SharedRAbundVectors& bv) : DataVector(bv), numBins(bv.numBins) {
        vector<RAbundVector*> data = bv.getSharedRAbundVectors();
        vector<string> temp = bv.getNamesGroups();
        for (int i = 0; i < temp.size(); i++) {
            push_back(data[i],temp[i]);
        }
        eliminateZeroOTUS();
        setLabel(bv.getLabel());
    }
    ~SharedRAbundVectors() { for (it = lookup.begin(); it != lookup.end(); it++) {  if (it->second != NULL) { delete it->second;  it->second = NULL; } }  lookup.clear(); }
    
    void setLabel(string l);
    int getOTUTotal(int bin);
    int get(int bin, string group);
    void set(int bin, int binSize, string group);
    int push_back(RAbundVector*, string group);
    int eliminateZeroOTUS(); //run after push_backs if groups are chosen
    int removeOTU(int bin);
    int removeOTUs(vector<int> bins);
    void removeGroups(vector<string> g);
    int size() { return lookup.size(); }
    int getNumBins() { return numBins; }
    int getNumSeqsSmallestGroup();
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedRAbundVectors.\n"); }
    void clear() { for (it = lookup.begin(); it != lookup.end(); it++) {  if (it->second != NULL) { delete it->second;  it->second = NULL; } }  lookup.clear(); }
    
    void print(ostream&);
    void printHeaders(ostream&);
    int getNumSeqs(string); //group
    
    RAbundVector getRAbundVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*);
    SharedOrderVector getSharedOrderVector();
    vector<string> getNamesGroups(); //same order as Rabunds
    vector<RAbundVector*> getSharedRAbundVectors();
    vector<RAbundFloatVector*> getSharedRAbundFloatVectors();
    
private:
    map<string, RAbundVector*> lookup;
    map<string, RAbundVector*>::iterator it;
    int numBins;
    
};


#endif /* sharedrabundvectors_hpp */
