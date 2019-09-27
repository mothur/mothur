//
//  accuracy.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "accuracy.hpp"

/***********************************************************************/
double Accuracy::getValue(double tp,  double tn,  double fp,  double fn) {
    try {
        long long p = tp + fn;
        long long n = fp + tn;
        double accuracy = (tp + tn) / (double) (p + n);
        if(p + n == 0)		{	accuracy = 0;								}
        
        if (isnan(accuracy) || isinf(accuracy)) { accuracy = 0; }
        
        return accuracy;
    }
    catch(exception& e) {
        m->errorOut(e, "Accuracy", "getValue");
        exit(1);
    }
}
/***********************************************************************/

