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

class SharedRAbundFloatVectors : public DataVector {
    
public:
    SharedRAbundFloatVectors() : DataVector() { label = ""; numBins = 0; }
    SharedRAbundFloatVectors(ifstream&);
    SharedRAbundFloatVectors(SharedRAbundFloatVectors& bv) : DataVector(bv), numBins(bv.numBins) {
        vector<RAbundFloatVector*> data = bv.getSharedRAbundFloatVectors();
        vector<string> temp = bv.getNamesGroups();
        for (int i = 0; i < temp.size(); i++) {
            push_back(data[i],temp[i]);
        }
        eliminateZeroOTUS();
        setLabel(bv.getLabel());
    }
    ~SharedRAbundFloatVectors() { for (it = lookup.begin(); it != lookup.end(); it++) {  if (it->second != NULL) { delete it->second;  it->second = NULL; } }  lookup.clear(); }
    
    float getOTUTotal(int bin);
    float removeOTUs(vector<int> bins);
    float removeOTU(int bin);
    void setLabel(string l);
    float get(int bin, string group);
    void set(int bin, float binSize, string group);
    int push_back(RAbundFloatVector*, string n);
    void print(ostream&);
    void printHeaders(ostream&);
    void removeGroups(vector<string> g);
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedRAbundVectors.\n"); }
    void clear() { for (it = lookup.begin(); it != lookup.end(); it++) {  if (it->second != NULL) { delete it->second;  it->second = NULL; } }  lookup.clear(); }
    int size() { return lookup.size();  }
    int getNumBins() { return numBins; }
    float getNumSeqs(string); //group
    float getNumSeqsSmallestGroup();
    
    vector<RAbundVector*> getSharedRAbundVectors();
    vector<RAbundFloatVector*> getSharedRAbundFloatVectors();
    RAbundVector getRAbundVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*);
    vector<string> getNamesGroups();
    int eliminateZeroOTUS(); //run after push_backs if groups are chosen
    
private:
    map<string, RAbundFloatVector*> lookup;
    map<string, RAbundFloatVector*>::iterator it;
    int numBins;

};

#endif /* sharedrabundfloatvectors_hpp */


