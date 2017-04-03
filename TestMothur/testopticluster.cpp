//
//  testopticluster.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/15/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "testopticluster.h"
#include "optimatrix.h"
#include "distancecommand.h"
#include "dataset.h"

/**************************************************************************************************/
TestOptiCluster::TestOptiCluster()  {  //setup
    m = MothurOut::getInstance();
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(100); //Fasta, name, group returned
    
    string inputString = "fasta=" + filenames[0];
    m->mothurOut("/******************************************/"); m->mothurOutEndLine();
    m->mothurOut("Running command: dist.seqs(" + inputString + ")"); m->mothurOutEndLine();
    m->mothurCalling = true;
    
    Command* distCommand = new DistanceCommand(inputString);
    distCommand->execute();
    
    map<string, vector<string> > outputFilenames = distCommand->getOutputFiles();
    
    delete distCommand;
    m->mothurCalling = false;
    
    columnFile = outputFilenames["column"][0];
    
    OptiMatrix testMatrix(columnFile, filenames[1], "name", "column", 0.03, false);
    
    setVariables(&testMatrix, "mcc");
}
/**************************************************************************************************/
TestOptiCluster::~TestOptiCluster() {
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); } //teardown
    m->mothurRemove(columnFile);
}
/**************************************************************************************************/
TEST(TestOptiCluster, myInitialize) {
    TestOptiCluster test;
    double initialMetricValue;
    
    EXPECT_EQ(0,(test.initialize(initialMetricValue, true, "singleton")));
}

TEST(TestOptiCluster, calcs) {
    TestOptiCluster test; long long tp,tn,fp,fn;
    
    tp=5000; tn=10000; fp=10; fn=200;
    ASSERT_NEAR(0.9694, test.calcMCC(tp,tn,fp,fn), 0.0001); //metric value
    
    tp=5000; tn=10000; fp=10; fn=200;
    ASSERT_NEAR(0.9615, test.calcSens(tp,tn,fp,fn), 0.0001);
    
    tp=5000; tn=10000; fp=10; fn=200;
    ASSERT_NEAR(0.9990, test.calcSpec(tp,tn,fp,fn), 0.0001);
    
    tp=5000; tn=10000; fp=10; fn=200;
    ASSERT_NEAR(0.9861, test.calcTPTN(tp,tn,fp,fn), 0.0001);
    
    tp=5000; tn=10000; fp=10; fn=200;
    ASSERT_NEAR(0.9861, test.calcFPFN(tp,tn,fp,fn), 0.0001);
}

TEST(TestOptiCluster, myMoveAdjustTFValues) {
    TestOptiCluster test;
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    double initialMetricValue;
    test.initialize(initialMetricValue, false, "singleton"); //no randomization
    
    ASSERT_NEAR(0.9700, test.moveAdjustTFValues(0, 10, 1, tp,tn,fp,fn), 0.0001); //metric value
}

TEST(TestOptiCluster, myUpdate) {
    TestOptiCluster test;
    double initialMetricValue;
    test.initialize(initialMetricValue, false, "singleton"); //no randomization
    test.update(initialMetricValue);
    
    ASSERT_NEAR(0.7853, initialMetricValue, 0.0001); //metric value
}
/**************************************************************************************************/
