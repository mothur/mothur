//
//  fpfn.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "fpfn.hpp"

/***********************************************************************/
double FPFN::getValue( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        long long p = fp + fn;
        
        double fpfn = 1.0 - (p / (double)(tp + tn + fp + fn)); //minimize
        
        if (isnan(fpfn) || isinf(fpfn)) { fpfn = 0; }
        
        return fpfn;
    }
    catch(exception& e) {
        m->errorOut(e, "FPFN", "getValue");
        exit(1);
    }
}
/***********************************************************************/

