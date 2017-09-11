//
//  sharedrabundfloatvector.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/25/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "sharedrabundfloatvector.hpp"

/***********************************************************************/

SharedRAbundFloatVector::SharedRAbundFloatVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0), group("") {}

/***********************************************************************/

SharedRAbundFloatVector::SharedRAbundFloatVector(int n) : DataVector(), data(n,0) , maxRank(0), numBins(n), numSeqs(0), group("") {}

/***********************************************************************/

SharedRAbundFloatVector::SharedRAbundFloatVector(vector<float> rav) :  DataVector(), maxRank(0), numBins(rav.size()), numSeqs(0), group("")  {
    try {
        data.assign(numBins, 0);
        for(int i=0;i<rav.size();i++){ set(i, rav[i]); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "SharedRAbundFloatVector");
        exit(1);
    }
}

/***********************************************************************/

SharedRAbundFloatVector::SharedRAbundFloatVector(vector<float> rav, int mr, int nb, int ns) :  DataVector(), group(""){
    try {
        numBins = nb;
        maxRank = mr;
        numSeqs = ns;
        data = rav;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "SharedRAbundFloatVector");
        exit(1);
    }
}
/***********************************************************************/


SharedRAbundFloatVector::SharedRAbundFloatVector(ifstream& f) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
    try {
        f >> label >> group >> numBins;
        
        data.assign(numBins, 0);
        float inputData;
        
        for(int i=0;i<numBins;i++){
            f >> inputData;
            set(i, inputData);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "SharedRAbundFloatVector");
        exit(1);
    }
}

/***********************************************************************/
SharedRAbundFloatVector::SharedRAbundFloatVector(ifstream& f, string l, string g, int n) : DataVector(), maxRank(0), numBins(n), numSeqs(0) {
    try {
        label = l;
        group = g;
        data.assign(numBins, 0);
        
        float inputData;
        for(int i=0;i<numBins;i++){
            f >> inputData;
            set(i, inputData);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "SharedRAbundFloatVector");
        exit(1);
    }
}

/***********************************************************************/

SharedRAbundFloatVector::~SharedRAbundFloatVector() { }

/***********************************************************************/

void SharedRAbundFloatVector::set(int binNumber, float newBinSize){
    try {
        int oldBinSize = data[binNumber];
        data[binNumber] = newBinSize;
        
        if(newBinSize > maxRank)	{	maxRank = newBinSize;	}
        
        numSeqs += (newBinSize - oldBinSize);
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "set");
        exit(1);
    }
}

/***********************************************************************/

float SharedRAbundFloatVector::get(int index){ return data[index]; }
/***********************************************************************/

void SharedRAbundFloatVector::clear(){
    numBins = 0;
    maxRank = 0;
    numSeqs = 0;
    group = "";
    data.clear();
}
/***********************************************************************/

void SharedRAbundFloatVector::push_back(float binSize){
    try {
        data.push_back(binSize);
        numBins++;
        
        if(binSize > maxRank){ maxRank = binSize; }
        
        numSeqs += binSize;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "push_back");
        exit(1);
    }
}
/***********************************************************************/

float SharedRAbundFloatVector::remove(int bin){
    try {
        float abund = data[bin];
        data.erase(data.begin()+bin);
        numBins--;
        
        if(abund == maxRank){ maxRank = m->max(data); }
        
        numSeqs -= abund;
        
        return abund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "remove");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVector::resize(int size){ data.resize(size); }
/***********************************************************************/
int SharedRAbundFloatVector::size(){ return data.size(); }
/***********************************************************************/
void SharedRAbundFloatVector::print(ostream& output){
    try {
        output << label;
        output << '\t' << group << '\t' << numBins;
        
        for(int i=0;i<numBins;i++){		output  << '\t' << data[i];		}
        output << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "nonSortedPrint");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundFloatVector::getNumBins(){
    return numBins;
}

/***********************************************************************/

float SharedRAbundFloatVector::getNumSeqs(){
    return numSeqs;
}

/***********************************************************************/

float SharedRAbundFloatVector::getMaxRank(){
    return maxRank;
}
/***********************************************************************/

RAbundVector SharedRAbundFloatVector::getRAbundVector(){
    RAbundVector rav;
    for(int i = 0; i < data.size(); i++) { rav.push_back(int(data[i])); }
    rav.setLabel(label);
    return rav;
}
/***********************************************************************/

SharedRAbundVector SharedRAbundFloatVector::getSharedRAbundVector(){
    SharedRAbundVector rav; rav.setLabel(label); rav.setGroup(group);
    for(int i = 0; i < data.size(); i++) { rav.push_back(int(data[i])); }
    return rav;
}
/***********************************************************************/

SAbundVector SharedRAbundFloatVector::getSAbundVector() {
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
        m->errorOut(e, "SharedRAbundFloatVector", "getSAbundVector");
        exit(1);
    }
}

/***********************************************************************/

OrderVector SharedRAbundFloatVector::getOrderVector(map<string,int>* nameMap = NULL) {
    try {
 m->mothurOut("[ERROR]: can not convert SharedRAbundVectors to an ordervector, ordervectors assume no zero OTUS.\n"); m->setControl_pressed(true);
 OrderVector o; return o; 
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVector", "getOrderVector");
        exit(1);
    }
}

/***********************************************************************/

