//
//  testsharedrabundvectors.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testsharedrabundvectors_hpp
#define testsharedrabundvectors_hpp

#include "sharedrabundvectors.hpp"
#include "gtest/gtest.h"


class TestSharedRabundVectors : public SharedRAbundVectors {
    
public:
    
    TestSharedRabundVectors();
    ~TestSharedRabundVectors();
    
    MothurOut* m;
    string sharedFile;
};



#endif /* testsharedrabundvectors_hpp */
