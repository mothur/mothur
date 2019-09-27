//
//  sensspeccalc.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/22/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef sensspeccalc_hpp
#define sensspeccalc_hpp

#include "mothurout.h"
#include "optimatrix.h"

class SensSpecCalc {
    
public:
    SensSpecCalc(OptiData& matrix, ListVector* list);
    ~SensSpecCalc(){}
    
    void getResults(OptiData& matrix, double& tp, double& tn, double& fp, double& fn);
    
private:
    Utils util;
    MothurOut* m;
    vector<vector< int> > otus;
};



#endif /* sensspeccalc_hpp */
