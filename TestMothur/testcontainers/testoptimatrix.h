//
//  testoptimatrix.h
//  Mothur
//
//  Created by Sarah Westcott on 6/6/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testoptimatrix__
#define __Mothur__testoptimatrix__

#include "optimatrix.h"

class TestOptiMatrix : public OptiMatrix {
    
public:
    
    TestOptiMatrix();
    ~TestOptiMatrix();
    
    using OptiMatrix::getCloseSeqs;
    using OptiMatrix::findDistFormat;
    using OptiMatrix::readPhylip;
    using OptiMatrix::readColumn;
    
    string columnFile;
    vector<string> filenames;
    
};

#endif /* defined(__Mothur__testoptimatrix__) */
