//
//  mcc.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "mcc.hpp"

/***********************************************************************/
double MCC::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
    
        double p = tp + fn;
        double n = fp + tn;
        double pPrime = tp + fp;
        double nPrime = tn + fn;
        
        double matthewsCorrCoef = ((tp * tn) - (fp * fn)) / sqrt(p * n * pPrime * nPrime);
        
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

