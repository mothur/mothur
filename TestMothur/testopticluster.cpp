//
//  testopticluster.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/15/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "testopticluster.h"
#include "catch.hpp"
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
}
/**************************************************************************************************/
TestOptiCluster::~TestOptiCluster() {
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); } //teardown
    m->mothurRemove(columnFile);
}
/**************************************************************************************************/
TEST_CASE("Testing OptiCluster Class") {
    TestOptiCluster testOcluster;
    OptiMatrix matrix(testOcluster.columnFile, testOcluster.filenames[1], "name", 0.03, false);
    testOcluster.setVariables(&matrix, "mcc");
    
    SECTION("Testing Initialize") {
        INFO("Using First 100 sequences of final.fasta and final.names") // Only appears on a FAIL
        
        double initialMetricValue;
        
        CAPTURE(testOcluster.initialize(initialMetricValue, true)); //
        
        CHECK(testOcluster.initialize(initialMetricValue, true) == 0); //metric value
    }
    
    SECTION("Testing calcMCC") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        
        CAPTURE(testOcluster.calcMCC(tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.calcMCC(tp,tn,fp,fn)*10000) == 9694); //metric value
    }
    
    SECTION("Testing calcSens") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        
        CAPTURE(testOcluster.calcSens(tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.calcSens(tp,tn,fp,fn)*10000) == 9615); //metric value
    }
    
    SECTION("Testing calcSpec") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        
        CAPTURE(testOcluster.calcSpec(tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.calcSpec(tp,tn,fp,fn)*10000) == 9990); //metric value
    }
    
    SECTION("Testing calcTPTN") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        
        CAPTURE(testOcluster.calcTPTN(tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.calcTPTN(tp,tn,fp,fn)*10000) == 9861); //metric value
    }
    
    SECTION("Testing calcTP2TN") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        
        CAPTURE(testOcluster.calcTP2TN(tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.calcTP2TN(tp,tn,fp,fn)*10000) == 16436); //metric value
    }
    
    SECTION("Testing calcFPFN") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        
        CAPTURE(testOcluster.calcFPFN(tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.calcFPFN(tp,tn,fp,fn)*10000) == 9861); //metric value
    }
    
    SECTION("Testing moveAdjustTFValues") {
        INFO("Using tp=5000, tn=10000, fp=10, fn=200 and mcc") // Only appears on a FAIL
        
        double tp,tn,fp,fn; tp=5000; tn=10000; fp=10; fn=200;
        double initialMetricValue;
        testOcluster.initialize(initialMetricValue, false); //no randomization
        
        CAPTURE(testOcluster.moveAdjustTFValues(0, 10, 1, tp,tn,fp,fn)); // Displays this variable on a FAIL
        
        CHECK((int)(testOcluster.moveAdjustTFValues(0, 10, 1, tp,tn,fp,fn)*10000) == 9700); //metric value
    }
    
    SECTION("Testing update") {
        INFO("Using mcc") // Only appears on a FAIL
        
        double initialMetricValue;
        testOcluster.initialize(initialMetricValue, false); //no randomization
        testOcluster.update(initialMetricValue);
        
        CAPTURE(initialMetricValue); // Displays this variable on a FAIL
        
        CHECK((int)(initialMetricValue*10000) == 7853); //metric value
    }
    
}
/**************************************************************************************************/
