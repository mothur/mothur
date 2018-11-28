//
//  testsubsample.cpp
//  Mothur
//
//  Created by Sarah Westcott on 11/15/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testsubsample.hpp"

/**************************************************************************************************/
TestSubSample::TestSubSample()  {  //setup
    m = MothurOut::getInstance();
}
/**************************************************************************************************/
TestSubSample::~TestSubSample() { }
/**************************************************************************************************/
TEST(Test_SubSample, getWeightedSample) {
    TestSubSample test;
    
    map<long long, long long> weights;
    weights[1] = 1;
    weights[2] = 5;
    weights[3] = 10;
    weights[4] = 15;
    weights[5] = 20;
    weights[6] = 25;
    weights[7] = 30;
    weights[8] = 35;
    weights[9] = 40;
    weights[10] = 45; //226 total reads represented
    
    set<long long> names = test.getWeightedSample(weights, 10); //select all the reads
    
    EXPECT_EQ(1,*names.begin());
    
    names = test.getWeightedSample(weights, 5); //select 1 the read
    
    set<long long>::iterator it = names.find(10);
    
    EXPECT_EQ((it != names.end()),true);
}

/**************************************************************************************************/
