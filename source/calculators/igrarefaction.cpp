//
//  igrarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/6/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "igrarefaction.hpp"

/***********************************************************************/
int compare_doubles(const void* a, const void* b)
{
    double* arg1 = (double *) a;
    double* arg2 = (double *) b;
    if( *arg1 < *arg2 ) return -1;
    else if( *arg1 == *arg2 ) return 0;
    else return 1;
}
/***********************************************************************/
double IGRarefaction::calcMu(t_IGParams *ptIGParams)
{
    double dMu = 0.0;
    
    ptIGParams->dC = coverage;
    
    DiversityUtils dutils("igrarefaction");
    
    dutils.solveF(1.0, 1.0e10, ptIGParams, 1.0e-7, &dMu);
    
    return dMu;
}
/***********************************************************************/
vector<double> IGRarefaction::getValues(SAbundVector* rank, vector<mcmcSample>& sampling){
    try {
        t_Data   tData;
        vector<double> results;
        
#ifdef USE_GSL
        
        DiversityUtils dutils("igrarefaction");
        
        dutils.loadAbundance(&tData, rank);
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl
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
        }
        
        adMu = (double *) malloc(sizeof(double)*nSamples);
        
        //printf("numSample = %d ",nSamples);
        for(int i = 0; i < nSamples; i++){
            adMu[i] = ((double) sampled)*calcMu(&atIGParams[i]);
            //printf("%f\n", adMu[i]);
            //fflush(stdout);
        }
        
        qsort(adMu, nSamples, sizeof(double), compare_doubles);
        
        dLower  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.025);
        dMedian = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.5);
        dUpper  = gsl_stats_quantile_from_sorted_data(adMu, 1, nSamples, 0.975);
        
        m->mothurOut("\nIGRarefaction - d_Lower = " + toString(dLower) + " d_Median = " + toString(dMedian) + " d_Upper = " + toString(dUpper) + "\n\n");
        
        if (isnan(dLower) || isinf(dLower))     { dLower = 0;   }
        if (isnan(dMedian) || isinf(dMedian))   { dMedian = 0;  }
        if (isnan(dUpper) || isinf(dUpper))     { dUpper = 0;   }
        
        results.push_back(dLower); results.push_back(dMedian); results.push_back(dUpper);
        
        free(adMu);
        
        dutils.freeAbundance(&tData);
#endif
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "IGRarefaction", "getValues");
        exit(1);
    }
}
/***********************************************************************/
