//
//  testopticluster.h
//  Mothur
//
//  Created by Sarah Westcott on 6/15/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__testopticluster__
#define __Mothur__testopticluster__


#include "gtest.h"
#include "opticluster.h"
#include "fakeoptimatrix.hpp"


class TestOptiCluster : public OptiCluster   {
    
public:
    
    TestOptiCluster();
    ~TestOptiCluster();
    
protected:
    MothurOut* m;
    ClusterMetric* metric;
    FakeOptiMatrix testMatrix;
    
    using OptiCluster::setVariables;
    using OptiCluster::initialize;
    using OptiCluster::update;
    using OptiCluster::getCloseFarCounts;
    
    FRIEND_TEST(TestOptiCluster, myInitialize);
    FRIEND_TEST(TestOptiCluster, getCloseFarCounts);
    FRIEND_TEST(TestOptiCluster, myUpdate);

};

#endif /* defined(__Mothur__testopticluster__) */
