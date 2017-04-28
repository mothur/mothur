//
//  testclustercalcs.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/18/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef testclustercalcs_hpp
#define testclustercalcs_hpp

#include "gtest.h"
#include "mcc.hpp"
#include "fakemcc.hpp"
#include "optimatrix.h"


class TestMCCCalc  : public MCC  {
    
public:
    
    TestMCCCalc();
    ~TestMCCCalc();
    
    MothurOut* m;
    FakeClusterCalcValues fake;
};


#endif /* testclustercalcs_hpp */
