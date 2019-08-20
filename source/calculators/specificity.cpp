//
//  specificity.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "specificity.hpp"

/***********************************************************************/
double Specificity::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        long long n = fp + tn;
        double specificity = tn / (double) n;
        
        if(n == 0)			{	specificity = 0;	}
        if (isnan(specificity) || isinf(specificity)) { specificity = 0; }
        
        return specificity;
    }
    catch(exception& e) {
        m->errorOut(e, "Specificity", "getValue");
        exit(1);
    }
}
/***********************************************************************/

