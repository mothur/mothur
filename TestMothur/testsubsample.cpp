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
    
    map<string, long long> weights;
    weights["seq1"] = 1;
    weights["seq2"] = 5;
    weights["seq3"] = 10;
    weights["seq4"] = 15;
    weights["seq5"] = 20;
    weights["seq6"] = 25;
    weights["seq7"] = 30;
    weights["seq8"] = 35;
    weights["seq9"] = 40;
    weights["seq10"] = 45; //226 total reads represented
    
    set<string> names = test.getWeightedSample(weights, 10); //select all the reads
    
    EXPECT_EQ("seq1",*names.begin());
    
    names = test.getWeightedSample(weights, 1); //select 1 the read
    
    //EXPECT_EQ("seq10",*names.begin());
}

/**************************************************************************************************/
