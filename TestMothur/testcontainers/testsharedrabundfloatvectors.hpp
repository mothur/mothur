//
//  testsharedrabundfloatvectors.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/14/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testsharedrabundfloatvectors_hpp
#define testsharedrabundfloatvectors_hpp

#include "sharedrabundfloatvectors.hpp"
#include "gtest/gtest.h"


class TestSharedRabundFloatVectors : public SharedRAbundFloatVectors {
    
public:
    
    TestSharedRabundFloatVectors();
    ~TestSharedRabundFloatVectors();
    
    MothurOut* m;
    string relabundFile;
};



#endif /* testsharedrabundfloatvectors_hpp */
