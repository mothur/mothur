//
//  testopticluster.h
//  Mothur
//
//  Created by Sarah Westcott on 6/15/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testopticluster__
#define __Mothur__testopticluster__

#include "opticluster.h"
#include "gtest/gtest.h"


class TestOptiCluster : public OptiCluster, ::testing::Test  {
    
public:
    
    TestOptiCluster();
    ~TestOptiCluster();
    
protected:
    MothurOut* m;
    string columnFile, phylipFile;
    vector<string> filenames;
    OptiMatrix* matrix;
    
    using OptiCluster::calcMCC;
    using OptiCluster::calcSens;
    using OptiCluster::calcSpec;
    using OptiCluster::calcTPTN;
    using OptiCluster::calcFPFN;
    using OptiCluster::moveAdjustTFValues;
    using OptiCluster::setVariables;
    using OptiCluster::initialize;
    using OptiCluster::update;
    
};

#endif /* defined(__Mothur__testopticluster__) */
