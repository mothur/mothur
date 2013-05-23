
#ifndef PDSSparCC_runSparcc_h
#define PDSSparCC_runSparcc_h

//
//  runSparcc.h
//  PDSSparCC
//
//  Created by Patrick Schloss on 10/31/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

/**************************************************************************************************/

//#include "sparcc.h"
#include "randomnumber.h"
#include "mothurout.h"

/**************************************************************************************************/

class CalcSparcc {
	
public:
	CalcSparcc(vector<vector<float> >, int, int, string);
    vector<vector<float> > getRho()    {   return median;  }
private:
    MothurOut* m;
    void addPseudoCount(vector<vector<float> >&);
    vector<float> getLogFractions(vector<vector<float> >, string);
    void getT_Matrix(vector<float>);
    
    
    void getT_Vector();
    void getD_Matrix();
    vector<float> getBasisVariances();
    vector<vector<float> > getBasisCorrelations(vector<float>);
    float getExcludedPairs(vector<vector<float> >, int&, int&);
    void excludeValues(int, int);
    void getMedian(vector<vector<vector<float> > >);

    vector<float> tMatrix;

    vector<vector<float> > dMatrix;
    vector<float> tVector;
    vector<vector<int> > excluded;
    vector<vector<float> > median;
    
    int numOTUs;
    int numGroups;
    string normalizationMethod;
    
    RandomNumberGenerator RNG;    
};

#endif

/**************************************************************************************************/
