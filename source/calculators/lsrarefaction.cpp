//
//  lsrarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/20/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "lsrarefaction.hpp"


/***********************************************************************/
LSRarefaction::LSRarefaction(double c) : coverage(c), DiversityCalculator(true) {}
/***********************************************************************/
int compare_doubles2(const void* a, const void* b)
{
    double* arg1 = (double *) a;
    double* arg2 = (double *) b;
    if( *arg1 < *arg2 ) return -1;
    else if( *arg1 == *arg2 ) return 0;
    else return 1;
}
/***********************************************************************/
vector<double> LSRarefaction::getValues(int numSeqs, vector<mcmcSample>& sampling){ //rank->getNumSeqs(); //nj
    try {
        
#ifdef USE_GSL
        
        DiversityUtils dutils("lsrarefaction");
        
        int nSamples = sampling.size();
        double*     adMu = nullptr;
        double dLower = 0.0, dMedian = 0.0, dUpper = 0.0;
        
        gsl_set_error_handler_off();
        
        t_LSParams* atLSParams;
        atLSParams = (t_LSParams *) malloc(nSamples*sizeof(t_LSParams)); //MAX_SAMPLES
        
        //load sampling data
        for (int i = 0; i < nSamples; i++) {
            if (m->getControl_pressed()) { free(atLSParams); return results; }
            
            atLSParams[i].dMDash = sampling[i].alpha;
            atLSParams[i].dV = sampling[i].beta;
            atLSParams[i].dNu = sampling[i].dNu;
            atLSParams[i].dC = coverage;
            atLSParams[i].n = 0;
        }
        
        adMu = (double *) malloc(sizeof(double)*nSamples);
        
        for(int i = 0; i < nSamples; i++){ adMu[i] = ((double) numSeqs)*dutils.calcMu(&atLSParams[i]); }
        
        qsort(adMu, nSamples, sizeof(double), compare_doubles2);
        
        dLower  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.025);
        dMedian = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.5);
        dUpper  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.975);
        
        if (isnan(dLower) || isinf(dLower))     { dLower = 0;   }
        if (isnan(dMedian) || isinf(dMedian))   { dMedian = 0;  }
        if (isnan(dUpper) || isinf(dUpper))     { dUpper = 0;   }
        
        m->mothurOut("\nLSRarefaction - d_Lower = " + toString(dLower) + " d_Median = " + toString(dMedian) + " d_Upper = " + toString(dUpper) + "\n\n");
        
        results.push_back(dLower); results.push_back(dMedian); results.push_back(dUpper);
        
        //printf("%.2e:%.2e:%.2e ", dLower, dMedian, dUpper);
        
        free(adMu);
        free(atLSParams);
        
#endif
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "LSRarefaction", "getValues");
        exit(1);
    }
}
/***********************************************************************/


