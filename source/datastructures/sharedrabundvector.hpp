//
//  sharedrabundvector.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/24/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef sharedrabundvector_hpp
#define sharedrabundvector_hpp

#include "datavector.hpp"
#include "ordervector.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"


/*  Data Structure for a rabund file.
 This class is a child to datavector.  It represents OTU information at a certain distance.
	A rabundvector can be converted into and ordervector, listvector or sabundvector.
	Each member of the internal container "data" represents an individual OTU.
	So data[0] = 6, because there are six member in that OTU.
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m
 rabundvector	=	6				3			2		1		1
 sabundvector	=	2		1		1		0		0		1
 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5
 
 a,b,g,j,l = sample A
 c,d,h,i   = sample B
 e,f,k,m   = sample C
 
 The sharedRabund class is very similar to rabund.  SharedRabund allows for 0 otus in a sample, rabunds do not. SharedRabund also know their group. SharedRabunds are stored as floats, but printed as integers or floats depending on whether it represents a shared or relabund file.
 sharedrabund   =   A 5 2 1 1 1 0
                    B 5 2 2 0 0 0
                    C 5 2 0 1 0 1
 */


//class RAbundFloatVector;
class OrderVector;

class SharedRAbundVector : public DataVector {
    
public:
    SharedRAbundVector();
    SharedRAbundVector(int);
    SharedRAbundVector(vector<int>, int, int, int);
    SharedRAbundVector(vector<int>);
    SharedRAbundVector(const SharedRAbundVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs), group(bv.group) { };
    SharedRAbundVector(ifstream&);
    SharedRAbundVector(ifstream& f, string l, string g, int); //filehandle, label, numBins
    ~SharedRAbundVector();
    
    int getNumBins();
    int getNumSeqs();
    int getMaxRank();
    
    int remove(int);
    void set(int, int);
    int get(int);
    vector<int> get() { return data; }
    int increment(int); //add 1 to bin
    void push_back(int);
    void resize(int);
    int size();
    void clear();
    
    void print(ostream&); //nonsorted
    
    RAbundVector getRAbundVector();
    RAbundFloatVector getRAbundFloatVector();
    SAbundVector getSAbundVector();
    OrderVector getOrderVector(map<string,int>*); 
    
    string getGroup() { return group; } //group = "" for rabunds without groupInfo
    void setGroup(string g) { group = g;  }
    
private:
    vector<int> data;
    int maxRank;
    int numBins;
    int numSeqs;
    string group;
};



#endif /* sharedrabundvector_hpp */
