//
//  igabundance.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "igabundance.hpp"

/***********************************************************************/

double IGAbundance::getValues(SAbundVector* rank) {
    try {
        double results = 0;
        
        if (isnan(results) || isinf(results)) { results= 0; }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "IGAbundance", "getValues");
        exit(1);
    }
}
/***********************************************************************/


