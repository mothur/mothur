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

class TestOptiCluster : public OptiCluster {
    
public:
    
    TestOptiCluster();
    ~TestOptiCluster();
    
    MothurOut* m;
    string columnFile, phylipFile;
    vector<string> filenames;
    
    using OptiCluster::calcMCC;
    using OptiCluster::calcSens;
    using OptiCluster::calcSpec;
    using OptiCluster::calcTPTN;
    using OptiCluster::calcTP2TN;
    using OptiCluster::calcFPFN;
    using OptiCluster::moveAdjustTFValues;
    using OptiCluster::setVariables;
    
};

#endif /* defined(__Mothur__testopticluster__) */
