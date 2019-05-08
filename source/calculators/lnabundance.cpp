//
//  lnabundace.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "lnabundance.hpp"

/***********************************************************************/

vector<double> LNAbundance::getValues(SAbundVector* rank, vector<mcmcSample>& sampling) {
    try {
        
        int maxRank = rank->getMaxRank(); //nMax
        maxRank = floor(pow(2.0,ceil(log((double) maxRank)/log(2.0)) + 2.0) + 1.0e-7);
        
        vector<double> results; results.resize(maxRank, 0.0);
        int nSamples = sampling.size();
        
        if (nSamples == 0) {  return results; }
        
#ifdef USE_GSL
        
        DiversityUtils dutils("lnabund");
        
        for(int i = 0; i < sampling.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            for (int j = 1; j <= maxRank; j++) {
                int nA = j;
                double dLog = 0.0, dP = 0.0;
                
                if(nA < 100){ //MAX_QUAD
                    dLog = dutils.logLikelihoodQuad(nA, sampling[i].alpha, sampling[i].beta); //nA, dMDash, dV
                }
                else{
                    dLog = dutils.logLikelihoodRampal(nA, sampling[i].alpha, sampling[i].beta); //nA, dMDash, dV
                }
                
                
                dP = exp(dLog);
                
                results[j - 1] += dP*sampling[i].ns;
            }
            
        }
        
        for (int i = 1; i<=maxRank; i++) {
            results[i-1] /= (double)nSamples;
            
            if (isnan(results[i-1]) || isinf(results[i-1])) { results[i-1] = 0.0; }
        }
        
#endif
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "LNAbundance", "getValues");
        exit(1);
    }
}
/***********************************************************************/


