//
//  tn.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "tn.hpp"

/***********************************************************************/
double TN::getValue( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        double tnmax = tn / (double)(tp + tn + fp + fn);
        
        if (isnan(tnmax) || isinf(tnmax)) { tnmax = 0; }
        
        return tnmax;
    }
    catch(exception& e) {
        m->errorOut(e, "TN", "getValue");
        exit(1);
    }
}
/***********************************************************************/

