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
TestClusterCalcs::TestClusterCalcs(string metricName) {  //setup
    if (metricName == "mcc")             { metric = new MCC();              }
    else if (metricName == "sens")       { metric = new Sensitivity();      }
    else if (metricName == "spec")       { metric = new Specificity();      }
    else if (metricName == "tptn")       { metric = new TPTN();             }
    else if (metricName == "tp")         { metric = new TP();               }
    else if (metricName == "tn")         { metric = new TN();               }
    else if (metricName == "fp")         { metric = new FP();               }
    else if (metricName == "fn")         { metric = new FN();               }
    else if (metricName == "f1score")    { metric = new F1Score();          }
    else if (metricName == "accuracy")   { metric = new Accuracy();         }
    else if (metricName == "ppv")        { metric = new PPV();              }
    else if (metricName == "npv")        { metric = new NPV();              }
    else if (metricName == "fdr")        { metric = new FDR();              }
    else if (metricName == "fpfn")       { metric = new FPFN();             }

}
/**************************************************************************************************/
TestClusterCalcs::~TestClusterCalcs() { delete metric; }
/**************************************************************************************************/

TEST(Test_Calc_ClusterCalcs, mcc) {
    TestClusterCalcs test("mcc");
    double result = test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn);
    ASSERT_NEAR(0.791646, result, 0.0001); //metric value
}

TEST(Test_Calc_ClusterCalcs, sens) {
    TestClusterCalcs test("sens");
    ASSERT_NEAR(0.699235, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, spec) {
    TestClusterCalcs test("spec");
    ASSERT_NEAR(0.999951, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, tptn) {
    TestClusterCalcs test("tptn");
    ASSERT_NEAR(0.9997691, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, tp) {
    TestClusterCalcs test("tp");
    ASSERT_NEAR(0.000423, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, tn) {
    TestClusterCalcs test("tn");
    ASSERT_NEAR(0.9993461, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, fp) {
    TestClusterCalcs test("fp");
    ASSERT_NEAR(0.999951, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, fn) {
    TestClusterCalcs test("fn");
    ASSERT_NEAR(0.999818, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, f1score) {
    TestClusterCalcs test("f1score");
    ASSERT_NEAR(0.7856801, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, accuracy) {
    TestClusterCalcs test("accuracy");
    ASSERT_NEAR(0.999769, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, ppv) {
    TestClusterCalcs test("ppv");
    ASSERT_NEAR(0.896514, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, npv) {
    TestClusterCalcs test("npv");
    ASSERT_NEAR(0.9998179, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, fdr) {
    TestClusterCalcs test("fdr");
    ASSERT_NEAR(0.896514, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

TEST(Test_Calc_ClusterCalcs, fpfn) {
    TestClusterCalcs test("fpfn");
    ASSERT_NEAR(0.999769, test.metric->getValue(test.fake.tp,test.fake.tn,test.fake.fp,test.fake.fn), 0.0001); //metric value
    
}

/**************************************************************************************************/
