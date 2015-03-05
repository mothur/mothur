//
//  compare.h
//  Mothur
//
//  Created by Sarah Westcott on 3/4/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__compare__
#define __Mothur__compare__

#include "mothurout.h"


class Compare {
public:
    
    Compare(){
        AA=0; AT=0; AG=0; AC=0;
        TA=0; TT=0; TG=0; TC=0;
        GA=0; GT=0; GG=0; GC=0;
        CA=0; CT=0; CG=0; CC=0;
        NA=0; NT=0; NG=0; NC=0;
        Ai=0; Ti=0; Gi=0; Ci=0; Ni=0;
        dA=0; dT=0; dG=0; dC=0;
        refName = "";
        queryName = "";
        weight = 1;
        matches = 0;
        mismatches = 0;
        total = 0;
        errorRate = 1.0000;
        sequence = "";
    }
    ~Compare(){}
    
    int AA, AT, AG, AC,	TA, TT, TG, TC,	GA, GT, GG, GC,	CA, CT, CG, CC,	NA, NT, NG, NC, Ai, Ti, Gi, Ci, Ni, dA, dT, dG, dC;
    string refName, queryName, sequence;
    double errorRate;
    int weight, matches, mismatches, total;

};


#endif /* defined(__Mothur__compare__) */
