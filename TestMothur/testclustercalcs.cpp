//
//  testclustercalcs.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/18/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "testclustercalcs.hpp"
#include "mcc.hpp"

/**************************************************************************************************/
TestMCCCalc::TestMCCCalc() : MCC() {  //setup
    m = MothurOut::getInstance();
}
/**************************************************************************************************/
TestMCCCalc::~TestMCCCalc() {}
/**************************************************************************************************/

TEST(TestClusterCalcs, mcc) {
    TestMCCCalc test;
    ASSERT_NEAR(0.791646, test.getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}
/**************************************************************************************************/
