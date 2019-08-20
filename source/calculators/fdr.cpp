//
//  fdr.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "fdr.hpp"

/***********************************************************************/
double FDR::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        long long pPrime = tp + fp;
        double falseDiscoveryRate = fp / (double) pPrime;
        
        if(pPrime == 0)		{	falseDiscoveryRate = 0;		}
        
        if (isnan(falseDiscoveryRate) || isinf(falseDiscoveryRate)) { falseDiscoveryRate = 0; }
        
        return (1.0-falseDiscoveryRate);

    }
    catch(exception& e) {
        m->errorOut(e, "FDR", "getValue");
        exit(1);
    }
}
/***********************************************************************/

