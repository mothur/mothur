//
//  testopticluster.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/15/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "testopticluster.h"
#include "distancecommand.h"
#include "dataset.h"

/**************************************************************************************************/
TestOptiCluster::TestOptiCluster()  {  //setup
    m = MothurOut::getInstance();
    metric = new MCC();
    setVariables(&testMatrix, metric);
}
/**************************************************************************************************/
TestOptiCluster::~TestOptiCluster() { delete metric; }
/**************************************************************************************************/
TEST(TestOptiCluster, myInitialize) {
    TestOptiCluster test;
    double initialMetricValue;
    
    EXPECT_EQ(0,(test.initialize(initialMetricValue, true, "singleton")));
}

TEST(TestOptiCluster, myUpdate) {
    TestOptiCluster test;
    double initialMetricValue;
    test.initialize(initialMetricValue, false, "singleton"); //no randomization
    test.update(initialMetricValue);
    
    //first round
    ASSERT_NEAR(1, initialMetricValue, 0.00001); //metric value
    
    test.update(initialMetricValue);
    
    //first round
    ASSERT_NEAR(1, initialMetricValue, 0.00001); //metric value
}

TEST(TestOptiCluster, getCloseFarCounts) {
    TestOptiCluster test;
    double initialMetricValue;
    test.initialize(initialMetricValue, false, "singleton"); //no randomization
    test.update(initialMetricValue);
    
    vector<long long> results = test.getCloseFarCounts(0, 31);
    
    ASSERT_EQ(results[0], 0); //number of close sequences in bin 31 to seq 0
    ASSERT_EQ(results[1], 10); //number of far sequences in bin 31 to seq 0
}
/**************************************************************************************************/
