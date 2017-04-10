//
//  fp.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "fp.hpp"

/***********************************************************************/
double FP::getValue( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        double fpmin = fp / (double)(tp + tn + fp + fn);
        
        if (isnan(fpmin) || isinf(fpmin)) { fpmin = 0; }
        
        return (1.0 - fpmin);
    }
    catch(exception& e) {
        m->errorOut(e, "FP", "getValue");
        exit(1);
    }
}
/***********************************************************************/

