//
//  lnrarefaction.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/13/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "lnrarefaction.hpp"

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
double f1_lnrarefaction(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
    double dTemp = (x - dMDash);
    double dExp  = x*((double) n) - exp(x) - 0.5*((dTemp*dTemp)/dV);
    double dRet  = exp(dExp);
    
    return dRet;
}
#ifdef USE_GSL
/***********************************************************************/
double logLikelihoodQuad(int n, double dMDash, double dV)
{
    gsl_integration_workspace *ptGSLWS =
    gsl_integration_workspace_alloc(1000);
    double dLogFac1   = 0.0, dLogFacN  = 0.0;
    double dResult = 0.0, dError = 0.0, dPrecision = 0.0;
    gsl_function tGSLF;
    t_LNParams tLNParams;
    double dEst = dMDash + ((double)n)*dV, dA = 0.0, dB = 0.0;
    
    tLNParams.n = n; tLNParams.dMDash = dMDash; tLNParams.dV = dV;
    
    tGSLF.function = &f1_lnrarefaction;
    tGSLF.params   = (void *) &tLNParams;
    
    dLogFac1 = log(2.0*M_PI*dV);
    
    if(n < 50){
        dLogFacN = gsl_sf_fact(n);
        dLogFacN = log(dLogFacN);
    }
    else{
        dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
    }
    
    DiversityUtils dutils("lnrarefaction");
    
    if(dEst > dV){
        double dMax = 0.0;
        double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
        double dVar   = 0.0;
        
        dutils.solveF(0.0, dUpper, (void *) &tLNParams, 1.0e-5, &dMax);
        
        dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
        
        dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
    }
    else{
        double dMax = 0.0;
        double dLower = dEst - dV;
        double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
        double dVar   = 0.0;
        
        dutils.solveF(dLower, dUpper, (void *) &tLNParams, 1.0e-5, &dMax);
        
        dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
        
        dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
    }
    
    if(n < 10){
        dPrecision = HI_PRECISION;
    }
    else{
        dPrecision = LO_PRECISION;
    }
    
    gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
    
    gsl_integration_workspace_free(ptGSLWS);
    
    return log(dResult) - dLogFacN -0.5*dLogFac1;
}
/***********************************************************************/
double fMu_LNRarefaction(double x, void* pvParams)
{
    t_IGParams* ptIGParams = (t_IGParams *) pvParams;
    
    double dMDD = ptIGParams->dAlpha + x;
    double dLogP0 = logLikelihoodQuad(0, dMDD, ptIGParams->dBeta);
    
    return (1.0 - exp(dLogP0)) - ptIGParams->dC;
}
/***********************************************************************/
int solveF_lnrarefaction(double x_lo, double x_hi, double (*f)(double, void*),
           void* params, double tol, double *xsolve)
{
    int status, iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    gsl_function F;
    
    F.function = f;
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
double calcMu_lnrarefaction(t_IGParams *ptLNParams)
{
    double dLogMu = 0.0;
    
    solveF_lnrarefaction(0, 1.0e7, fMu_LNRarefaction, ptLNParams, 1.0e-7, &dLogMu);
    
    return exp(dLogMu);
}
#endif
/***********************************************************************/
vector<double> LNRarefaction::getValues(SAbundVector* rank, vector<mcmcSample>& sampling){
    try {
        t_Data   tData;
        vector<double> results;
        
#ifdef USE_GSL
        
        DiversityUtils dutils("lnrarefaction");
        
        dutils.loadAbundance(&tData, rank);
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl
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
        
        for(int i = 0; i < nSamples; i++){ adMu[i] = ((double) sampled)*calcMu_lnrarefaction(&atIGParams[i]); }
        
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
        
        dutils.freeAbundance(&tData);
#endif
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "LNRarefaction", "getValues");
        exit(1);
    }
}
/***********************************************************************/
