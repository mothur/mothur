//
//  lnshift.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/14/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "lnshift.hpp"

/***********************************************************************/
LNShift::LNShift() : DiversityCalculator(true) {}
/***********************************************************************/

vector<double> LNShift::getValues(int numSeqs, vector<mcmcSample>& sampling) {
    try {
        
        int nMax = 100000; //nMax
        
        results.resize(nMax, 0.0);
        int nSamples = sampling.size();
        
        if (nSamples == 0) {  return results; }
        
#ifdef USE_GSL
        
        DiversityUtils dutils("lnshift");
        
        gsl_set_error_handler_off();
        
        double dShift = log(5.0e5/(double)numSeqs);
        
        for(int i = 0; i < sampling.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            for (int j = 1; j <= nMax; j++) {
                int nA = j;
                double dLog = 0.0, dP = 0.0;
                
                if(nA < 100){ //MAX_QUAD
                    dLog = dutils.logLikelihoodQuad(nA, sampling[i].alpha + dShift, sampling[i].beta);
                }
                else{
                    dLog = dutils.logLikelihoodRampal(nA, sampling[i].alpha + dShift, sampling[i].beta);
                }
                
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
        m->errorOut(e, "LNShift", "getValues");
        exit(1);
    }
}
/***********************************************************************/

