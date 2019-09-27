//
//  ppv.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "ppv.hpp"

/***********************************************************************/
double PPV::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        long long pPrime = tp + fp;
        double positivePredictiveValue = tp / (double) pPrime;
        
        if(pPrime == 0)		{	positivePredictiveValue = 0;		}
        
        if (isnan(positivePredictiveValue) || isinf(positivePredictiveValue)) { positivePredictiveValue = 0; }
        
        return positivePredictiveValue;
    }
    catch(exception& e) {
        m->errorOut(e, "PPV", "getValue");
        exit(1);
    }
}
/***********************************************************************/

