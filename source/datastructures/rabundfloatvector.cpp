//
//  rabundfloatvector.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "rabundfloatvector.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "ordervector.hpp"


/***********************************************************************/

RAbundFloatVector::RAbundFloatVector() : DataVector(), maxRank(0.0), numBins(0), numSeqs(0.0), group("") {}

/***********************************************************************/

RAbundFloatVector::RAbundFloatVector(int n) : DataVector(), data(n,0) , maxRank(0), numBins(0), numSeqs(0), group("") {}

/***********************************************************************/

//RAbundVector::RAbundVector(const RAbundVector& rav) : DataVector(rav), data(rav.data), (rav.label),  (rav.maxRank), (rav.numBins), (rav.numSeqs){}


/***********************************************************************/

RAbundFloatVector::RAbundFloatVector(string id, vector<float> rav) : DataVector(id), data(rav), group("") {
    try {
        numBins = 0;
        maxRank = 0;
        numSeqs = 0;
        
        for(int i=0;i<data.size();i++){
            if(!util.isEqual(data[i], 0))		{	numBins = i+1;		}
            if(data[i] > maxRank)	{	maxRank = data[i];	}
            numSeqs += data[i];
        }
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "RAbundFloatVector");
        exit(1);
    }
}

/***********************************************************************/

RAbundFloatVector::RAbundFloatVector(vector<float> rav, float mr, int nb, float ns) : group(""){
    try {
        numBins = nb;
        maxRank = mr;
        numSeqs = ns;
        data = rav;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "RAbundFloatVector");
        exit(1);
    }
}
/***********************************************************************/


RAbundFloatVector::RAbundFloatVector(ifstream& f) : DataVector(), maxRank(0), numBins(0), numSeqs(0), group("") {
    try {
        int hold;
        f >> label >> hold;
        
        data.assign(hold, 0);
        float inputData;
        
        for(int i=0;i<hold;i++){
            f >> inputData;
            set(i, inputData);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "RAbundFloatVector");
        exit(1);
    }
}

/***********************************************************************/
RAbundFloatVector::RAbundFloatVector(ifstream& f, string l, string g) : DataVector(), maxRank(0), numBins(0), numSeqs(0), group(g) {
    try {
        int hold;
        label = l;
        f >> hold;
        
        float inputData;
        for(int i=0;i<hold;i++){
            f >> inputData;
            push_back(inputData);
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "RAbundFloatVector");
        exit(1);
    }
}

/***********************************************************************/

RAbundFloatVector::~RAbundFloatVector() {}

/***********************************************************************/

void RAbundFloatVector::set(int binNumber, float newBinSize){
    try {
        float oldBinSize = data[binNumber];
        data[binNumber] = newBinSize;
        
        if(util.isEqual(oldBinSize, 0))			{	numBins++;				}
        if(util.isEqual(newBinSize, 0))			{	numBins--;				}
        if(newBinSize > maxRank)	{	maxRank = newBinSize;	}
        
        numSeqs += (newBinSize - oldBinSize);
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "set");
        exit(1);
    }
}

/***********************************************************************/

float RAbundFloatVector::get(int index){
    return data[index];
    
}
/***********************************************************************/

void RAbundFloatVector::clear(){
    numBins = 0;
    maxRank = 0;
    numSeqs = 0;
    data.clear();
    
}
/***********************************************************************/

void RAbundFloatVector::push_back(float binSize){
    try {
        data.push_back(binSize);
        numBins++;
        
        if(binSize > maxRank){
            maxRank = binSize;
        }
        
        numSeqs += binSize;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "push_back");
        exit(1);
    }
}

/***********************************************************************/

void RAbundFloatVector::pop_back(){
    
    return data.pop_back();
}

/***********************************************************************/

void RAbundFloatVector::resize(int size){
    
    data.resize(size);
}

/***********************************************************************/

int RAbundFloatVector::size(){
    return data.size();
}

/***********************************************************************/

void RAbundFloatVector::quicksort(){
    sort(data.rbegin(), data.rend());
}

/***********************************************************************/

float RAbundFloatVector::sum(){
    return sum(0);
}

/***********************************************************************/

float RAbundFloatVector::sum(int index){
    float sum = 0;
    for(int i = index; i < data.size(); i++) {  sum += data[i];  }
    return sum;
}
/***********************************************************************/

float RAbundFloatVector::remove(int bin){
    try {
        float abund = data[bin];
        data.erase(data.begin()+bin);
        numBins--;
        
        if(util.isEqual(abund, maxRank)){ maxRank = util.max(data); }
        
        numSeqs -= abund;
        
        return abund;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundVector", "remove");
        exit(1);
    }
}

/***********************************************************************/

int RAbundFloatVector::numNZ(){
    int numNZ = 0;
    for(int i = 0; i < data.size(); i++) { if(!util.isEqual(data[i], 0)) { numNZ++; } }
    return numNZ;
}

/***********************************************************************/

vector<float>::reverse_iterator RAbundFloatVector::rbegin(){
    return data.rbegin();
}

/***********************************************************************/

vector<float>::reverse_iterator RAbundFloatVector::rend(){
    return data.rend();
}

/***********************************************************************/
void RAbundFloatVector::nonSortedPrint(ostream& output){
    try {
        output << label;
        if (group != "") { output << '\t' << group; }
        output << '\t' << numBins;

        
        for(int i=0;i<numBins;i++){		output  << '\t' << data[i];		}
        output << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "nonSortedPrint");
        exit(1);
    }
}
/***********************************************************************/
void RAbundFloatVector::print(ostream& output){
    try {
        output << label;
        if (group != "") { output << '\t' << group; }
        output << '\t' << numBins;
        
        vector<float> hold = data;
        sort(hold.rbegin(), hold.rend());
        
        for(int i=0;i<numBins;i++){		output  << '\t' << hold[i];		}
        output << endl;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "print");
        exit(1);
    }
}
/***********************************************************************/
int RAbundFloatVector::getNumBins(){
    return numBins;
}

/***********************************************************************/

float RAbundFloatVector::getNumSeqs(){
    return numSeqs;
}

/***********************************************************************/

float RAbundFloatVector::getMaxRank(){
    return maxRank;
}

/***********************************************************************/

RAbundFloatVector RAbundFloatVector::getRAbundFloatVector(){
    return *this;
}
/***********************************************************************/

RAbundVector RAbundFloatVector::getRAbundVector(){
    RAbundVector rav;
    rav.setLabel(label);
    for (int i = 0; i < data.size(); i++) { rav.push_back((int)data[i]); }
    return rav;
}
/***********************************************************************/

SAbundVector RAbundFloatVector::getSAbundVector() {
    try {
        SAbundVector sav(maxRank+1);
        
        for(int i=0;i<data.size();i++){
            int abund = (int)data[i];
            sav.set(abund, sav.get(abund) + 1);
        }
        sav.set(0, 0);
        sav.setLabel(label);
        return sav;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "getSAbundVector");
        exit(1);
    }
}

/***********************************************************************/

OrderVector RAbundFloatVector::getOrderVector(map<string,int>* nameMap) {
    try {
        OrderVector ov;
        
        for(int i=0;i<data.size();i++){
            for(int j=0;j<data[i];j++){
                ov.push_back(i);
            }
        }
        util.mothurRandomShuffle(ov);
        ov.setLabel(label);	
        ov.getNumBins();
        
        return ov;
    }
    catch(exception& e) {
        m->errorOut(e, "RAbundFloatVector", "getOrderVector");
        exit(1);
    }
}

/***********************************************************************/
