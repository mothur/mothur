//
//  testsharedrabundvector.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/9/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testsharedrabundvector_hpp
#define testsharedrabundvector_hpp

#include "sharedrabundvector.hpp"
#include "gtest/gtest.h"


class TestSharedRabundVector : public SharedRAbundVector {
    
public:
    
    TestSharedRabundVector();
    ~TestSharedRabundVector();
    
    MothurOut* m;
    string sharedFile;
};


#endif /* testsharedrabundvector_hpp */
