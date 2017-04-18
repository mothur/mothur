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
#include "optimatrix.h"


class TestOptiCluster : public OptiCluster   {
    
public:
    
    TestOptiCluster();
    ~TestOptiCluster();
    
protected:
    MothurOut* m;
    string columnFile, phylipFile;
    vector<string> filenames;
    ClusterMetric* metric;
    OptiMatrix* testMatrix;
    
    using OptiCluster::setVariables;
    using OptiCluster::initialize;
    using OptiCluster::update;
    
    FRIEND_TEST(TestOptiCluster, myInitialize);
    FRIEND_TEST(TestOptiCluster, myMoveAdjustTFValues);
    FRIEND_TEST(TestOptiCluster, myUpdate);
    FRIEND_TEST(TestOptiCluster, calcs);

};

#endif /* defined(__Mothur__testopticluster__) */
