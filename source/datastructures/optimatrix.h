//
//  optimatrix.h
//  Mothur
//
//  Created by Sarah Westcott on 4/20/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__optimatrix__
#define __Mothur__optimatrix__

#include "optidata.hpp"
/*
#ifdef UNIT_TEST
    OptiMatrix() : OptiData(0.03) {};
#endif
*/

class OptiMatrix : public OptiData {

public:
    OptiMatrix(vector< set<long long> >, vector<string>, vector<string>, double); //closeness, namemap, singleton, cutoff
    OptiMatrix(string, string, string, string, double, bool); //distfile, dupsFile, dupsFormat, distFormat, cutoff, sim
    ~OptiMatrix(){}
    
protected:
    
    string distFile, namefile, countfile, format, distFormat;
    bool sim;

    int readPhylip();
    int readColumn();
};


#endif /* defined(__Mothur__optimatrix__) */
