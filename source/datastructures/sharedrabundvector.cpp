//
//  sharedrabundvector.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/24/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "sharedrabundvector.hpp"


/***********************************************************************/

SharedRAbundVector::SharedRAbundVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0), group("") {}

/***********************************************************************/

SharedRAbundVector::SharedRAbundVector(int n) : DataVector(), data(n,0) , maxRank(0), numBins(n), numSeqs(0), group("") {}

/***********************************************************************/

SharedRAbundVector::SharedRAbundVector(vector<int> rav) :  DataVector(), maxRank(0), numBins(rav.size()), numSeqs(0), group("")  {
    try {
        data.assign(numBins, 0);
        for(int i=0;i<rav.size();i++){ set(i, rav[i]); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "SharedRAbundVector");
        exit(1);
    }
}

/***********************************************************************/

SharedRAbundVector::SharedRAbundVector(vector<int> rav, int mr, int nb, int ns) :  DataVector(), group(""){
    try {
        numBins = nb;
        maxRank = mr;
        numSeqs = ns;
        data = rav;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "SharedRAbundVector");
        exit(1);
    }
}
/***********************************************************************/


SharedRAbundVector::SharedRAbundVector(ifstream& f) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
    try {
        f >> label >> group >> numBins;
        
        data.assign(numBins, 0);
        int inputData;
        
        for(int i=0;i<numBins;i++){
            f >> inputData;
            set(i, inputData);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "SharedRAbundVector");
        exit(1);
    }
}

/***********************************************************************/
SharedRAbundVector::SharedRAbundVector(ifstream& f, string l, string g, int n) : DataVector(), maxRank(0), numBins(n), numSeqs(0) {
    try {
        label = l;
        group = g;
        data.assign(numBins, 0);
        
        int inputData;
        for(int i=0;i<numBins;i++){
            f >> inputData;
            set(i, inputData);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "SharedRAbundVector");
        exit(1);
    }
}

/***********************************************************************/

SharedRAbundVector::~SharedRAbundVector() { }

/***********************************************************************/

void SharedRAbundVector::set(int binNumber, int newBinSize){
    try {
        int oldBinSize = data[binNumber];
        data[binNumber] = newBinSize;
        
        if(newBinSize > maxRank)	{	maxRank = newBinSize;	}
        
        numSeqs += (newBinSize - oldBinSize);
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundVector", "set");
        exit(1);
    }
}
/***********************************************************************/

int SharedRAbundVector::increment(int binNumber){
    try {
        data[binNumber]++;
        int newBinSize = data[binNumber];
        
        if(newBinSize > maxRank)	{	maxRank = newBinSize;	}
        
        numSeqs++;
        
        return newBinSize;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundVector", "increment");
        exit(1);
    }
}

/***********************************************************************/

int SharedRAbundVector::get(int index){ return data[index]; }
/***********************************************************************/

void SharedRAbundVector::clear(){
    numBins = 0;
    maxRank = 0;
    numSeqs = 0;
    group = "";
    data.clear();
}
/***********************************************************************/

void SharedRAbundVector::push_back(int binSize){
    try {
        data.push_back(binSize);
        numBins++;
        
        if(binSize > maxRank){ maxRank = binSize; }
        
        numSeqs += binSize;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "push_back");
        exit(1);
    }
}
/***********************************************************************/

int SharedRAbundVector::remove(int bin){
    try {
        int abund = data[bin];
        data.erase(data.begin()+bin);
        numBins--;
        
        if(abund == maxRank){
            vector<int>::iterator it = max_element(data.begin(), data.end());
            maxRank = *it;
        }
        
        numSeqs -= abund;
        
        return abund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "remove");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundVector::resize(int size){ data.resize(size); }
/***********************************************************************/
int SharedRAbundVector::size(){ return data.size(); }
/***********************************************************************/
void SharedRAbundVector::print(ostream& output){
    try {
        output << label;
        output << '\t' << group << '\t' << numBins;
        
        for(int i=0;i<numBins;i++){		output  << '\t' << int(data[i]);		}
        output << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "nonSortedPrint");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVector::getNumBins(){
    return numBins;
}

/***********************************************************************/

int SharedRAbundVector::getNumSeqs(){
    return numSeqs;
}

/***********************************************************************/

int SharedRAbundVector::getMaxRank(){
    return maxRank;
}
/***********************************************************************/

RAbundVector SharedRAbundVector::getRAbundVector(){
    RAbundVector rav;
    for(int i = 0; i < data.size(); i++) { if (data[i] != 0) {  rav.push_back(int(data[i])); } }
    rav.setLabel(label);
    return rav;
}
/***********************************************************************/

SAbundVector SharedRAbundVector::getSAbundVector() {
    try {
        SAbundVector sav(int(maxRank+1));
        
        for(int i=0;i<data.size();i++){
            int abund = data[i];
            sav.set(abund, int(sav.get(abund)) + 1);
        }
        sav.set(0, 0);
        sav.setLabel(label);
        return sav;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundVector", "getSAbundVector");
        exit(1);
    }
}

/***********************************************************************/

OrderVector SharedRAbundVector::getOrderVector(map<string,int>* nameMap = NULL) {
    try {
    m->mothurOut("[ERROR]: can not convert SharedRAbundVectors to an ordervector, ordervectors assume no zero OTUS.\n"); m->setControl_pressed(true);
            OrderVector o; return o;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundVector", "getOrderVector");
        exit(1);
    }
}

/***********************************************************************/
