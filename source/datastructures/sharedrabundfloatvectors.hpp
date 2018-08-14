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
#include "sharedrabundfloatvector.hpp"
#include "sharedordervector.h"
#include "ordervector.hpp"

/*  DataStructure for a relabund file. */
//********************************************************************************************************************
inline bool compareRAbundFloats(SharedRAbundFloatVector* left, SharedRAbundFloatVector* right){ return (left->getGroup() < right->getGroup()); }
//********************************************************************************************************************


class SharedRAbundFloatVectors : public DataVector {
    
public:
    SharedRAbundFloatVectors() : DataVector() { label = ""; numBins = 0;  otuTag = "Otu"; printSharedHeaders = true; }
    SharedRAbundFloatVectors(ifstream&, vector<string>&, string&, string&);
    SharedRAbundFloatVectors(SharedRAbundFloatVectors& bv) : DataVector(bv), numBins(bv.numBins), otuTag(bv.otuTag) {
        vector<SharedRAbundFloatVector*> data = bv.getSharedRAbundFloatVectors();
        for (int i = 0; i < data.size(); i++) { push_back(data[i]); }
        setLabels(bv.getLabel());
        printSharedHeaders = true;
        setOTUNames(bv.getOTUNames());
        eliminateZeroOTUS();
    }
    ~SharedRAbundFloatVectors() { clear(); }
    
    void setLabels(string l);
    float getOTUTotal(int bin);
    vector<float> getOTU(int bin);
    float removeOTU(int bin);
    float get(int bin, string group);
    void set(int bin, float binSize, string group);
    void setOTUNames(vector<string> names);
    vector<string> getOTUNames();
    string getOTUName(int);
    void setOTUName(int, string);

    int push_back(SharedRAbundFloatVector*);
    void removeGroups(vector<string> g);
    int removeGroups(int minSize, bool silent=false);  // removes any groups with numSeqs < minSize
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedRAbundFloatVectors.\n"); }
    void clear() { for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear(); groupNames.clear(); numBins = 0; currentLabels.clear(); }
    int size() { return lookup.size();  }
    int getNumGroups() { return lookup.size(); }
    int getNumBins() { return numBins; }
    float getNumSeqs(string); //group
    float getNumSeqsSmallestGroup();
    void print(ostream&);
    
    vector<SharedRAbundVector*> getSharedRAbundVectors();
    vector<SharedRAbundFloatVector*> getSharedRAbundFloatVectors();
    RAbundVector getRAbundVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*) { m->mothurOut("[ERROR]: can not convert SharedRAbundVectors to an ordervector, ordervectors assume no zero OTUS.\n"); m->setControl_pressed(true); OrderVector o; return o; }
    
    vector<string> getNamesGroups();
    void eliminateZeroOTUS(); //run after push_backs if groups are chosen
    
private:
    void printHeaders(ostream&);
    vector<SharedRAbundFloatVector*> lookup;
    vector<string> currentLabels;
    map<string, int> groupNames;
    int numBins;
    string otuTag;
    bool printSharedHeaders;

};

#endif /* sharedrabundfloatvectors_hpp */


