//
//  tptn.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "tptn.hpp"

/***********************************************************************/
double TPTN::getValue(double tp,  double tn,  double fp,  double fn)  {
    try {
        long long p = tp + tn;
        double tptn = p / (double)(tp + tn + fp + fn);
        
        if (isnan(tptn) || isinf(tptn)) { tptn = 0; }
        
        return tptn;
    }
    catch(exception& e) {
        m->errorOut(e, "TPTN", "getValue");
        exit(1);
    }
}
/***********************************************************************/

