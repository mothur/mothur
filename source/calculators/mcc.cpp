//
//  mcc.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "mcc.hpp"

/***********************************************************************/
double MCC::getValue( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        
        long long p = tp + fn;
        long long n = fp + tn;
        long long pPrime = tp + fp;
        long long nPrime = tn + fn;
        
        double matthewsCorrCoef = ((tp * tn) - (fp * fn)) / (double) sqrt(p * n * pPrime * nPrime);
        if(p == 0 || n == 0 || pPrime == 0 || nPrime == 0){	matthewsCorrCoef = 0;	}
        
        if (isnan(matthewsCorrCoef) || isinf(matthewsCorrCoef)) { matthewsCorrCoef = 0; }
        
        return matthewsCorrCoef;
    }
    catch(exception& e) {
        m->errorOut(e, "MCC", "getValue");
        exit(1);
    }
}
/***********************************************************************/

