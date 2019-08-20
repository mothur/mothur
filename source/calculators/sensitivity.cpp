//
//  sensitivity.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "sensitivity.hpp"

/***********************************************************************/
double Sensitivity::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        
        long long p = tp + fn;
        double sensitivity = tp / (double) p;
        
        if(p == 0)	{	sensitivity = 0;	}
        if (isnan(sensitivity) || isinf(sensitivity)) { sensitivity = 0; }
        
        return sensitivity;
    }
    catch(exception& e) {
        m->errorOut(e, "Sensitivity", "getValue");
        exit(1);
    }
}
/***********************************************************************/

