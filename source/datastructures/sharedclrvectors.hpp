//
//  sharedclrvectors.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/21/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef sharedclrvectors_hpp
#define sharedclrvectors_hpp

#include "datavector.hpp"
#include "sharedclrvector.hpp"

/*  DataStructure for a clr relabund file.
 
 The clr - log centered ratio - is the log (base 2) of a value divided by the geometric mean of the values across all OTUs for that sample. For example here are the counts for four OTUs in one sample...

 > x <- c(10, 5, 3, 1)
 > log2(x / prod(x)^(1/4))
 [1]  1.5147234  0.5147234 -0.2222422 -1.8072047
 
 */
//********************************************************************************************************************
inline bool compareCLRVectors(SharedCLRVector* left, SharedCLRVector* right){ return (left->getGroup() < right->getGroup()); }
//********************************************************************************************************************


class SharedCLRVectors : public DataVector {
    
public:
    SharedCLRVectors() : DataVector() { label = ""; numBins = 0;  otuTag = "Otu";  }
    SharedCLRVectors(ifstream&, vector<string>&, string&, string&);
    SharedCLRVectors(SharedCLRVectors& bv) : DataVector(bv), numBins(bv.numBins), otuTag(bv.otuTag) {
        vector<SharedCLRVector*> data = bv.getSharedCLRVectors();
        for (int i = 0; i < data.size(); i++) { push_back(data[i]); }
        setLabels(bv.getLabel());
        setOTUNames(bv.getOTUNames());
       //eliminateZeroOTUS();
    }
    ~SharedCLRVectors() { clear(); }
    
    vector<string> getNamesGroups();
    void setLabels(string l);
    float getOTUTotal(int bin);
    vector<float> getOTU(int bin);
    float removeOTU(int bin);
    float removeOTUs(vector<int>, bool sorted=false); //bins to remove, sorted or not
    float get(int bin, string group);
    void set(int bin, float binSize, string group);
    void setOTUNames(vector<string> names);
    vector<string> getOTUNames();
    string getOTUName(int);
    void setOTUName(int, string);

    int push_back(vector<float>, string binLabel=""); //add otu. mothur assumes abunds are in same order as groups.
    int push_back(SharedCLRVector*);
    void removeGroups(vector<string> g);
    int removeGroups(int minSize, bool silent=false);  // removes any groups with numSeqs < minSize
    void resize(int n) { m->mothurOut("[ERROR]: can not use resize for SharedCLRVectors.\n"); }
    void clear() { for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear(); groupNames.clear(); numBins = 0; currentLabels.clear(); }
    int size() { return (int)lookup.size();  }
    int getNumGroups() { return (int)lookup.size(); }
    int getNumBins() { return numBins; }
    float getNumSeqs(string); //group
    float getNumSeqsSmallestGroup();
    void print(ostream&, bool&);
    
    vector<SharedCLRVector*> getSharedCLRVectors();
    
    RAbundVector getRAbundVector() { m->mothurOut("[ERROR]: can not use getRAbundVector for SharedCLRVectors.\n"); RAbundVector r; return r; }
    SAbundVector getSAbundVector() { m->mothurOut("[ERROR]: can not use getSAbundVector for SharedCLRVectors.\n"); SAbundVector s; return s; }
    OrderVector getOrderVector(map<string,int>* hold = NULL) { m->mothurOut("[ERROR]: can not use getOrderVector for SharedCLRVectors.\n"); OrderVector o; return o; }
    
private:
    void printHeaders(ostream&, bool&);
    vector<SharedCLRVector*> lookup;
    vector<string> currentLabels;
    map<string, int> groupNames;
    int numBins;
    string otuTag;
};

#endif /* sharedclrvectors_hpp */
