//
//  igrarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/6/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "igrarefaction.hpp"

/***********************************************************************/
inline int compare_doubles1(const void* a, const void* b)
{
    double* arg1 = (double *) a;
    double* arg2 = (double *) b;
    if( *arg1 < *arg2 ) return -1;
    else if( *arg1 == *arg2 ) return 0;
    else return 1;
}
#ifdef USE_GSL
/***********************************************************************/
double fMu_igrarefaction(double x, void* pvParams)
{
    t_IGParams* ptIGParams = (t_IGParams *) pvParams;
    // double tx = x / 1667.0;
    
    DiversityUtils dutils("igrarefaction");
    
    double dAlphaDD = ptIGParams->dAlpha*sqrt(x);
    double dBetaDD  = ptIGParams->dBeta*x;
    double dLogP0   = dutils.logLikelihood(0, dAlphaDD, dBetaDD);
    
    // printf("dAlpha %f dBeta %f x %f dAlphaDD %f  dBetaDD %f dLofP0 %f", ptIGParams->dAlpha, ptIGParams->dBeta, x, dAlphaDD, dBetaDD, dLogP0);
    
    return (1.0 - exp(dLogP0)) - ptIGParams->dC;
}
/***********************************************************************/
int solveF_igrarefaction(double x_lo, double x_hi, void* params, double tol, double *xsolve)
{
    int status, iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    gsl_function F;
    
    F.function = &fMu_igrarefaction;
    F.params = params;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    
    do{
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        
        status = gsl_root_test_interval (x_lo, x_hi, 0, tol);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    (*xsolve) = gsl_root_fsolver_root (s);
    gsl_root_fsolver_free (s);
    
    return status;
}
/***********************************************************************/
 double calcMu_igrarefaction(t_IGParams *ptIGParams)
 {
     double dMu = 0.0;
 
     solveF_igrarefaction(1.0, 1.0e10, ptIGParams, 1.0e-7, &dMu);
  
     return dMu;
 }
#endif
/***********************************************************************/
vector<double> IGRarefaction::getValues(SAbundVector* rank, vector<mcmcSample>& sampling){
    try {
        t_Data   tData;
        vector<double> results;
        
#ifdef USE_GSL
        
        DiversityUtils dutils("igrarefaction");
        
        dutils.loadAbundance(&tData, rank);
        
        int sampled = rank->getNumSeqs(); //nj
        int nSamples = sampling.size();
        
        double*     adMu = NULL;
        double dLower = 0.0, dMedian = 0.0, dUpper = 0.0;
        
        gsl_set_error_handler_off();
        
        t_IGParams* atIGParams;
        atIGParams = (t_IGParams *) malloc(nSamples*sizeof(t_IGParams)); //MAX_SAMPLES
        atIGParams->dC = coverage;
        
        //load sampling data
        for (int i = 0; i < nSamples; i++) {
            if (m->getControl_pressed()) { return results; }
            
            atIGParams[i].dAlpha = sampling[i].alpha;
            atIGParams[i].dBeta = sampling[i].beta;
            atIGParams[i].nS = sampling[i].ns;
        }
        
        adMu = (double *) malloc(sizeof(double)*nSamples);
        
        for(int i = 0; i < nSamples; i++){ adMu[i] = ((double) sampled)*calcMu_igrarefaction(&atIGParams[i]); }
        
        qsort(adMu, nSamples, sizeof(double), compare_doubles1);
        
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
