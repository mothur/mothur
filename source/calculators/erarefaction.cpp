//
//  erarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "erarefaction.hpp"

/***********************************************************************/
ERarefaction::ERarefaction(int inc) : increment(inc), DiversityCalculator(true) {}
/***********************************************************************/

int ERarefaction::getValues(SAbundVector* rank, vector<double>& values){
    try {
        int maxRank = rank->getMaxRank();
        int sampled = rank->getNumSeqs(); //nl
        int numOTUs = rank->getNumBins(); //ns
        
        for (int n = 1; n <= sampled; n++) {
            if((n % increment) == 0){
                
                double dSum = 0.0;
#ifdef USE_GSL
                double dDenom = gsl_sf_lnchoose(sampled, n);
                
                for(int i = 1; i <= maxRank; i++){
                    
                    if (m->getControl_pressed()) { break; }
                    
                    int abund = rank->get(i);
                    
                    if (abund != 0) {
                        int thisRank = i; //nA
                        
                        if(sampled - thisRank >= n){
                            
                            double dNumer = gsl_sf_lnchoose(sampled - thisRank, n);
                            
                            dSum += ((double) abund)*exp(dNumer - dDenom);
                        }
                    }
                }
#endif
                double result = ((double) numOTUs) - dSum;
                
                if (isnan(result) || isinf(result)) { result = 0; }
                
                values.push_back(result);

            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ERarefaction", "getValues");
        exit(1);
    }
}
/***********************************************************************/
