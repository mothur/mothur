//
//  lnrarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/13/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "lnrarefaction.hpp"

/***********************************************************************/
LNRarefaction::LNRarefaction(double c) : coverage(c), DiversityCalculator(true) {}
/***********************************************************************/
inline int compare_doubles(const void* a, const void* b)
{
    double* arg1 = (double *) a;
    double* arg2 = (double *) b;
    if( *arg1 < *arg2 ) return -1;
    else if( *arg1 == *arg2 ) return 0;
    else return 1;
}
/***********************************************************************/
vector<double> LNRarefaction::getValues(int numSeqs, vector<mcmcSample>& sampling){ //rank->getNumSeqs(); //nj
    try {
        
#ifdef USE_GSL
        
        DiversityUtils dutils("lnrarefaction");
        
        int nSamples = sampling.size();
        double*     adMu = NULL;
        double dLower = 0.0, dMedian = 0.0, dUpper = 0.0;
        
        gsl_set_error_handler_off();
        
        t_IGParams* atIGParams;
        atIGParams = (t_IGParams *) malloc(nSamples*sizeof(t_IGParams)); //MAX_SAMPLES
        
        //load sampling data
        for (int i = 0; i < nSamples; i++) {
            if (m->getControl_pressed()) { return results; }
            
            atIGParams[i].dAlpha = sampling[i].alpha;
            atIGParams[i].dBeta = sampling[i].beta;
            atIGParams[i].nS = sampling[i].ns;
            atIGParams[i].dC = coverage;
        }
        
        adMu = (double *) malloc(sizeof(double)*nSamples);
        
        for(int i = 0; i < nSamples; i++){ adMu[i] = ((double) numSeqs)*dutils.calcMu(&atIGParams[i]); }
        
        qsort(adMu, nSamples, sizeof(double), compare_doubles);
        
        dLower  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.025);
        dMedian = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.5);
        dUpper  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.975);
        
        m->mothurOut("\nLNRarefaction - d_Lower = " + toString(dLower) + " d_Median = " + toString(dMedian) + " d_Upper = " + toString(dUpper) + "\n\n");
        
        if (isnan(dLower) || isinf(dLower))     { dLower = 0;   }
        if (isnan(dMedian) || isinf(dMedian))   { dMedian = 0;  }
        if (isnan(dUpper) || isinf(dUpper))     { dUpper = 0;   }
        
        results.push_back(dLower); results.push_back(dMedian); results.push_back(dUpper);
        
        free(adMu);
#endif
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "LNRarefaction", "getValues");
        exit(1);
    }
}
/***********************************************************************/
