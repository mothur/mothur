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


class OptiMatrix : public OptiData {
    
#ifdef UNIT_TEST
    friend class TestOptiMatrix;
    friend class FakeOptiMatrix;
#endif

    
public:
    
    OptiMatrix(vector< set<int> >, vector<string>, vector<string>, double);
    OptiMatrix(string, string, string, string, double, bool);
    ~OptiMatrix(){ }
    
protected:
    
    string distFile, namefile, countfile, format, distFormat;
    bool sim;

    int readPhylip();
    int readColumn();
};


#endif /* defined(__Mothur__optimatrix__) */
