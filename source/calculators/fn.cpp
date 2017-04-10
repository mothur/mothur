//
//  fn.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "fn.hpp"

/***********************************************************************/
double FN::getValue( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        double fnmin = fn / (double)(tp + tn + fp + fn);
        
        if (isnan(fnmin) || isinf(fnmin)) { fnmin = 0; }
        
        return (1.0 - fnmin);
    }
    catch(exception& e) {
        m->errorOut(e, "FN", "getValue");
        exit(1);
    }
}
/***********************************************************************/

