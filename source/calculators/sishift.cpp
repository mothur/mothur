//
//  sishift.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "sishift.hpp"

/***********************************************************************/

vector<double> SIShift::getValues(SAbundVector* rank, vector<mcmcSample>& sampling) {
    try {
        
        int nMax = 1000; //nMax
        nMax = floor(pow(2.0,ceil(log((double) nMax)/log(2.0)) + 2.0) + 1.0e-7);
        
        vector<double> results; results.resize(nMax, 0.0);
        int nSamples = (int)sampling.size();
        
        if (nSamples == 0) {  return results; }
        
#ifdef USE_GSL
        
        DiversityUtils dutils("sishift");
        
        gsl_set_error_handler_off();
        
        double dShift = 5.0e5/(double) rank->getNumSeqs();
        
        for(int i = 0; i < nSamples; i++) {
            
            if (m->getControl_pressed()) { break; }
            
            for (int j = 1; j <= nMax; j++) {
                int nA = j;
                double dLog = 0.0, dP = 0.0;
            
                dLog = dutils.logLikelihood(nA, sampling[i].alpha*sqrt(dShift), sampling[i].beta*dShift, sampling[i].dNu);
                
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
        m->errorOut(e, "SIShift", "getValues");
        exit(1);
    }
}
/***********************************************************************/

