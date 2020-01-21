//
//  sharedlcrvector.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/21/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef sharedlcrvector_hpp
#define sharedlcrvector_hpp

#include "datavector.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "ordervector.hpp"

class SharedLCRVector : public DataVector {
    
public:
    SharedLCRVector();
    SharedLCRVector(int);
    SharedLCRVector(vector<float>, float, int, float); //maxRank, numbins, numSeqs
    SharedLCRVector(vector<float>);
    SharedLCRVector(const SharedLCRVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs), group(bv.group) {};
    SharedLCRVector(ifstream&);
    SharedLCRVector(ifstream& f, string l, string g, int); //filehandle, label
    ~SharedLCRVector(){}
    
    int getNumBins()        {  return numBins;  }
    float getNumSeqs()      {  return numSeqs;  }
    float getMaxRank()      {  return maxRank;  }
    
    float remove(int);
    float remove(vector<int>);
    void set(int, float);
    float get(int);
    vector<float> get() { return data; }
    void push_back(float);
    void resize(int);
    int size();
    void clear();
    
    void print(ostream&); //nonsorted

    string getGroup() { return group; } //group = "" for rabunds without groupInfo
    void setGroup(string g) { group = g;  }
    
    RAbundVector getRAbundVector() { m->mothurOut("[ERROR]: can not use getRAbundVector for SharedLCRVector.\n"); RAbundVector r; return r; }
    SAbundVector getSAbundVector() { m->mothurOut("[ERROR]: can not use getSAbundVector for SharedLCRVector.\n"); SAbundVector s; return s; }
    OrderVector getOrderVector(map<string,int>* hold = NULL) { m->mothurOut("[ERROR]: can not use getOrderVector for SharedLCRVector.\n"); OrderVector o; return o; }
    
private:
    vector<float> data;
    float maxRank;
    int numBins;
    float numSeqs;
    string group;
};


#endif /* sharedlcrvector_hpp */
