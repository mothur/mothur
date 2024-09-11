//
//  testsharedrabundfloatvector.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/9/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testsharedrabundfloatvector_hpp
#define testsharedrabundfloatvector_hpp

#include "sharedrabundfloatvector.hpp"
#include "gtest.h"


class TestSharedRabundFloatVector : public SharedRAbundFloatVector {
    
public:
    
    TestSharedRabundFloatVector();
    ~TestSharedRabundFloatVector();
    
    MothurOut* m;
    string relabundFile;
};



#endif /* testsharedrabundfloatvector_hpp */
