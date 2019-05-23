//
//  siabundance.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/22/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "siabundance.hpp"

/***********************************************************************/
SIAbundance::SIAbundance() : DiversityCalculator(true) {}
/***********************************************************************/

vector<double> SIAbundance::getValues(int nMax, vector<mcmcSample>& sampling) { //int nMax = rank->getMaxRank();
    try {
        
        nMax = floor(pow(2.0,ceil(log((double) nMax)/log(2.0)) + 2.0) + 1.0e-7);
        
        vector<double> results; results.resize(nMax, 0.0);
        int nSamples = sampling.size();
        
        if (nSamples == 0) {  return results; }
        
#ifdef USE_GSL
        
        DiversityUtils dutils("siabund");
        
        gsl_set_error_handler_off();
        
        for(int i = 0; i < sampling.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            for (int j = 1; j <= nMax; j++) {
                int nA = j;
                double dLog = 0.0, dP = 0.0;
                
                dLog = dutils.logLikelihood(nA, sampling[i].alpha, sampling[i].beta, sampling[i].dNu);
               
                dP = exp(dLog);
                
                results[j - 1] += dP*sampling[i].ns;
            }
        }
        
        for (int i = 1; i<=nMax; i++) {
            results[i-1] /= (double)nSamples;
            
            if (isnan(results[i-1]) || isinf(results[i-1])) { results[i-1] = 0.0; }
        }
        
#endif
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "SIAbundance", "getValues");
        exit(1);
    }
}
/***********************************************************************/

