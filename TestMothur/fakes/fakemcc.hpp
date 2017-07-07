//
//  fakemcc.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/18/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fakemcc_hpp
#define fakemcc_hpp

#include "mothurout.h"


class FakeClusterCalcValues {
    
public:
    
    FakeClusterCalcValues() {  tp = 823; tn = 1944106; fp = 95; fn = 354; }
    ~FakeClusterCalcValues() {}
    
    long long tp, tn, fp, fn;
};



#endif /* fakemcc_hpp */

