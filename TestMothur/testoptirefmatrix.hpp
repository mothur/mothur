//
//  testoptirefmatrix.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/24/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testoptirefmatrix_hpp
#define testoptirefmatrix_hpp

#include "optirefmatrix.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"


class TestOptiRefMatrix : public OptiRefMatrix {
    
public:
    
    TestOptiRefMatrix();
    ~TestOptiRefMatrix();
    
    string columnFile, phylipFile;
    vector<string> filenames;
    vector<string> reffilenames;
    
private:
    
    
};


#endif /* testoptirefmatrix_hpp */
