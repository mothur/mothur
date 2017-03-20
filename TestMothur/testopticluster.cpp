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
TestOptiCluster::TestOptiCluster() {  //setup
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
    
    matrix = new OptiMatrix(columnFile, filenames[1], "name", "column", 0.03, false);
    
    setVariables(matrix, "mcc");
}
/**************************************************************************************************/
TestOptiCluster::~TestOptiCluster() {
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); } //teardown
    m->mothurRemove(columnFile);
    delete matrix;
}
/**************************************************************************************************/
TEST_F(TestOptiCluster, myInitialize) {
    double initialMetricValue;
    
    EXPECT_EQ(0,(initialize(initialMetricValue, true, "singleton")));
}

TEST_F(TestOptiCluster, calcMCC) {
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    
    ASSERT_DOUBLE_EQ(0.9694, calcMCC(tp,tn,fp,fn)); //metric value
}

TEST_F(TestOptiCluster, calcSens) {
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    
    ASSERT_DOUBLE_EQ(0.9615, calcSens(tp,tn,fp,fn)); //metric value
}

TEST_F(TestOptiCluster, calcSpec) {
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    
    ASSERT_DOUBLE_EQ(0.9990, calcSpec(tp,tn,fp,fn)); //metric value
}

TEST_F(TestOptiCluster, calcTPTN) {
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    
    ASSERT_DOUBLE_EQ(0.9861, calcTPTN(tp,tn,fp,fn)); //metric value
}

TEST_F(TestOptiCluster, calcFPFN) {
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    
    ASSERT_DOUBLE_EQ(0.9861, calcFPFN(tp,tn,fp,fn)); //metric value
}

TEST_F(TestOptiCluster, moveAdjustTFValues) {
    long long tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
    double initialMetricValue;
    initialize(initialMetricValue, false, "singleton"); //no randomization
    
    ASSERT_DOUBLE_EQ(0.9700, moveAdjustTFValues(0, 10, 1, tp,tn,fp,fn)); //metric value
}

TEST_F(TestOptiCluster, update) {
    double initialMetricValue;
    initialize(initialMetricValue, false, "singleton"); //no randomization
    update(initialMetricValue);
    
    ASSERT_DOUBLE_EQ(0.7853, initialMetricValue); //metric value
}
/**************************************************************************************************/
