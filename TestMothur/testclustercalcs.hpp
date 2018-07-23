//
//  testclustercalcs.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/18/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef testclustercalcs_hpp
#define testclustercalcs_hpp

#include "gtest/gtest.h"
#include "mcc.hpp"
#include "sensitivity.hpp"
#include "specificity.hpp"
#include "fdr.hpp"
#include "npv.hpp"
#include "ppv.hpp"
#include "f1score.hpp"
#include "tp.hpp"
#include "fp.hpp"
#include "fpfn.hpp"
#include "tptn.hpp"
#include "tn.hpp"
#include "fn.hpp"
#include "accuracy.hpp"
#include "fakemcc.hpp"


class TestClusterCalcs   {
    
public:
    
    TestClusterCalcs(string);
    ~TestClusterCalcs();
    
    FakeClusterCalcValues fake;
    ClusterMetric* metric;
    
private:
    MothurOut* m;
    
   
};


#endif /* testclustercalcs_hpp */
