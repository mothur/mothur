//
//  testsubsample.hpp
//  Mothur
//
//  Created by Sarah Westcott on 11/15/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testsubsample_hpp
#define testsubsample_hpp

#include "gtest/gtest.h"
#include "subsample.h"

class TestSubSample : public SubSample   {
    
public:
    
    TestSubSample();
    ~TestSubSample();
    
    MothurOut* m;
};

#endif /* testsubsample_hpp */
