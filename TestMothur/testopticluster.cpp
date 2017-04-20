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
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(); //Fasta, name, group returned
    columnFile = data.getSubsetFNGDistFile();
    
    testMatrix = new OptiMatrix(columnFile, filenames[1], "name", "column", 0.03, false);
    
    metric = new MCC();
    
    setVariables(testMatrix, metric);
}
/**************************************************************************************************/
TestOptiCluster::~TestOptiCluster() {
    delete metric; delete testMatrix;
}
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
    ASSERT_NEAR(0.95524, initialMetricValue, 0.00001); //metric value
}
/**************************************************************************************************/
