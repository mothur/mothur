//
//  tp.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "tp.hpp"

/***********************************************************************/
double TP::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        double tpmax = tp / (double)(tp + tn + fp + fn);
        
        if (isnan(tpmax) || isinf(tpmax)) { tpmax = 0; }
        
        return tpmax;
    }
    catch(exception& e) {
        m->errorOut(e, "TP", "getValue");
        exit(1);
    }
}
/***********************************************************************/

