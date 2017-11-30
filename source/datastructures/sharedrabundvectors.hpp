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
#include "sharedrabundvector.hpp"
#include "sharedrabundfloatvector.hpp"

/*  DataStructure for a shared file. */

//********************************************************************************************************************
inline bool compareRAbunds(SharedRAbundVector* left, SharedRAbundVector* right){ return (left->getGroup() < right->getGroup()); }
//********************************************************************************************************************

class SharedRAbundVectors : public DataVector {
    
public:
    SharedRAbundVectors() : DataVector() {  label = ""; numBins = 0; otuTag = "Otu"; printSharedHeaders = true; }
    SharedRAbundVectors(ifstream&, vector<string>& userGroups, string&, string&);
    SharedRAbundVectors(SharedRAbundVectors& bv) : DataVector(bv), numBins(bv.numBins), otuTag(bv.otuTag) {
        printSharedHeaders = true;
        vector<SharedRAbundVector*> data = bv.getSharedRAbundVectors();
        for (int i = 0; i < data.size(); i++) { push_back(data[i]); }
        eliminateZeroOTUS();
        setLabels(bv.getLabel());
    }
    ~SharedRAbundVectors() { clear(); }
    
    void setLabels(string l);
    int getOTUTotal(int bin);
    vector<int> getOTU(int bin);
    int removeOTU(int bin);
    int push_back(vector<int>, string binLabel=""); //add otu. mothur assumes abunds are in same order as groups.
    int get(int bin, string group);
    void set(int bin, int binSize, string group);
    int push_back(SharedRAbundVector*);
    void setOTUNames(vector<string> names);
    vector<string> getOTUNames();
    string getOTUName(int);
    void setOTUName(int, string);
    
    void eliminateZeroOTUS(); //run after push_backs if groups are chosen
    
    void removeGroups(vector<string> g);
    int removeGroups(int minSize, bool silent=false);  // removes any groups with numSeqs < minSize
    int size() { return lookup.size(); }
    int getNumBins() { return numBins; }
    int getNumSeqsSmallestGroup();
    vector<string> getNamesGroups(); //same order as Rabunds
    int getNumGroups() { return lookup.size(); }
    int getNumSeqs(string); //group
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedRAbundVectors.\n"); m->setControl_pressed(true); }
    void clear() { for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear(); groupNames.clear(); numBins = 0; }
    
    void print(ostream&);
    void printHeaders(ostream&);
    
    RAbundVector getRAbundVector();
    RAbundVector getRAbundVector(string); //group you want the rabund for
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*) { m->mothurOut("[ERROR]: can not convert SharedRAbundVectors to an ordervector, ordervectors assume no zero OTUS.\n"); m->setControl_pressed(true); OrderVector o; return o; }
    SharedOrderVector getSharedOrderVector();
    vector<SharedRAbundVector*> getSharedRAbundVectors();
    vector<SharedRAbundFloatVector*> getSharedRAbundFloatVectors();
    
private:
    vector<SharedRAbundVector*> lookup;
    vector<string> currentLabels;
    map<string, int> groupNames;
    int numBins;
    string otuTag;
    bool printSharedHeaders;
    
};


#endif /* sharedrabundvectors_hpp */
