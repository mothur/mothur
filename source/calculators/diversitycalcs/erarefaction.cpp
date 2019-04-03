//
//  erarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "erarefaction.hpp"

/***********************************************************************/

EstOutput ERarefaction::getValues(SAbundVector* rank){
    try {
        data = 0;
        
        int maxRank = rank->getMaxRank();
        int sampled = rank->getNumSeqs(); //nl
        int numOTUs = rank->getNumBins(); //ns
        
        double dDenom = gsl_sf_lnchoose(sampled, 1);
        double dSum = 0.0;
        
        for(int i = 1; i < maxRank; i++){
            
            if (m->getControl_pressed()) { break; }
            
            int abund = rank->get(i);
            if (abund != 0) {
                int thisRank = i; //nA
                if(sampled - thisRank >= 1){
                    
                    double dNumer = gsl_sf_lnchoose(sampled - thisRank, 1);
                    
                    dSum += ((double) abund)*exp(dNumer - dDenom);
                }
            }
        }
        
        data = ((double) numOTUs) - dSum;

        if (isnan(data) || isinf(data)) { data = 0; }
        
        return data;
    }
    catch(exception& e) {
        m->errorOut(e, "ERarefaction", "getValues");
        exit(1);
    }
}
/***********************************************************************/
