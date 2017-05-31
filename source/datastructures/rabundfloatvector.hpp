//
//  rabundfloatvector.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef rabundfloatvector_hpp
#define rabundfloatvector_hpp

#include "datavector.hpp"

/*  Data Structure for a rabund file.
 This class is a child to datavector.  It represents OTU information at a certain distance.
	A rabundvector can be converted into and ordervector, listvector or sabundvector.
	Each member of the internal container "data" represents an individual OTU.
	So data[0] = 6, because there are six member in that OTU.
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m
 rabundvector	=	6				3			2		1		1
 sabundvector	=	2		1		1		0		0		1
 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */

class RAbundFloatVector : public DataVector {
    
public:
    RAbundFloatVector();
    RAbundFloatVector(int);
    RAbundFloatVector(vector<float>, float, int, float);
    RAbundFloatVector(string, vector<float>);
    RAbundFloatVector(const RAbundFloatVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs){};
    RAbundFloatVector(ifstream&);
    RAbundFloatVector(ifstream& f, string l); //label given
    ~RAbundFloatVector();
    
    int getNumBins();
    float getNumSeqs();
    float getMaxRank();
    
    void set(int, float);
    float get(int);
    void push_back(float);
    float remove(int);
    void pop_back();
    void resize(int);
    int size();
    void quicksort();
    float sum();
    float sum(int);
    int numNZ();
    void clear();
    vector<float>::reverse_iterator rbegin();
    vector<float>::reverse_iterator rend();
    
    void print(ostream&); //sorted, no group
    int print(ostream&, string); //notsorted, group
    void print(string, ostream&); //label, sorted
    void nonSortedPrint(ostream&); //nonsorted , no group
    
    RAbundFloatVector getRAbundFloatVector();
    RAbundVector getRAbundVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*);
    
private:
    vector<float> data;
    float maxRank;
    int numBins;
    float numSeqs;
};



#endif /* rabundfloatvector_hpp */
