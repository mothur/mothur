//
//  npv.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "npv.hpp"

/***********************************************************************/
double NPV::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        long long nPrime = tn + fn;
        double negativePredictiveValue = tn / (double) nPrime;
        
        if(nPrime == 0)		{	negativePredictiveValue = 0;		}
        
        if (isnan(negativePredictiveValue) || isinf(negativePredictiveValue)) { negativePredictiveValue = 0; }
        
        return negativePredictiveValue;
    }
    catch(exception& e) {
        m->errorOut(e, "NPV", "getValue");
        exit(1);
    }
}
/***********************************************************************/

