//
//  testopticluster.h
//  Mothur
//
//  Created by Sarah Westcott on 6/15/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testopticluster__
#define __Mothur__testopticluster__


#include "gtest/gtest.h"
#include "opticluster.h"
#include "fakeoptimatrix.hpp"


class TestOptiCluster : public OptiCluster   {
    
public:
    
    TestOptiCluster();
    ~TestOptiCluster();
    

    MothurOut* m;
    ClusterMetric* metric;
    FakeOptiMatrix testMatrix;
    
    using OptiCluster::setVariables;
    using OptiCluster::initialize;
    using OptiCluster::update;
    using OptiCluster::getCloseFarCounts;
};

#endif /* defined(__Mothur__testopticluster__) */
