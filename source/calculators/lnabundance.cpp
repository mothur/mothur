//
//  lnabundace.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "lnabundance.hpp"

double f1_0(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
    double dTemp = (x - dMDash);
    double dExp  = x*((double) n) - exp(x) - 0.5*((dTemp*dTemp)/dV);
    double dRet  = exp(dExp);
    return dRet;
}

#ifdef USE_GSL
double logLikelihoodRampal2(int n, double dMDash, double dV)
{
    double dN = (double) n;
    double dLogLik = 0.0, dTemp = gsl_pow_int(log(dN) - dMDash,2), dTemp3 = gsl_pow_int(log(dN) - dMDash,3);
    
    dLogLik = -0.5*log(2.0*M_PI*dV) - log(dN) - (dTemp/(2.0*dV));
    
    dLogLik += log(1.0 + 1.0/(2.0*dN*dV)*(dTemp/dV + log(dN) - dMDash - 1.0)
                   + 1.0/(6.0*dN*dN*dV*dV*dV)*(3.0*dV*dV - (3.0*dV - 2.0*dV*dV)*(dMDash - log(dN))
                                               - 3.0*dV*dTemp + dTemp3));
    
    return dLogLik;
}
double derivExponent1(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
    double dTemp = (x - dMDash)/dV, dRet = 0.0;
    
    dRet = ((double) n) - exp(x) - dTemp;
    
    return dRet;
}

int solveF(double x_lo, double x_hi, double (*f)(double, void*),
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

double logLikelihoodQuad2(int n, double dMDash, double dV)
{
    gsl_integration_workspace *ptGSLWS =
    gsl_integration_workspace_alloc(1000);
    double dLogFac1   = 0.0, dLogFacN  = 0.0;
    double dResult = 0.0, dError = 0.0, dPrecision = 0.0;
    gsl_function tGSLF;
    t_LNParams tLNParams;
    double dEst = dMDash + ((double)n)*dV, dA = 0.0, dB = 0.0;
    
    tLNParams.n = n; tLNParams.dMDash = dMDash; tLNParams.dV = dV;
    
    tGSLF.function = &f1_0;
    tGSLF.params   = (void *) &tLNParams;
    
    dLogFac1 = log(2.0*M_PI*dV);
    
    
    //DiversityUtils dutils("lnabund");
    
    if(n < 50){
        dLogFacN = gsl_sf_fact(n);
        dLogFacN = log(dLogFacN);
    }
    else{
        dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
    }
    
    if(dEst > dV){
        double dMax = 0.0;
        double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
        double dVar   = 0.0;
        
        solveF(0.0, dUpper, derivExponent1, (void *) &tLNParams, 1.0e-5, &dMax);
        
        dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
        
        dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
    }
    else{
        double dMax = 0.0;
        double dLower = dEst - dV;
        double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
        double dVar   = 0.0;
        
        solveF(dLower, dUpper,derivExponent1,  (void *) &tLNParams, 1.0e-5, &dMax);
        
        dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
        
        dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
    }
    
    if(n < 10){
        dPrecision = HI_PRECISION;
    }
    else{
        dPrecision = LO_PRECISION;
    }
    printf("dA = %f dB = %f dPrecision = %f, GSL_INTEG_GAUSS61 %d maxlevel %zu \n", dA, dB, dPrecision, GSL_INTEG_GAUSS61, ptGSLWS->maximum_level);
    printf("alist = %f blist = %f  rlist = %f  elist = %f  order = %zu  level = %zu \n", *ptGSLWS->alist, *ptGSLWS->blist, *ptGSLWS->rlist, *ptGSLWS->elist, *ptGSLWS->order, *ptGSLWS->level);
    
    gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
    
    gsl_integration_workspace_free(ptGSLWS);
    
    return log(dResult) - dLogFacN -0.5*dLogFac1;
}
#endif
/***********************************************************************/

vector<double> LNAbundance::getValues(SAbundVector* rank, vector<mcmcSample>& sampling) {
    try {
        
        int maxRank = rank->getMaxRank(); //nMax
        maxRank = floor(pow(2.0,ceil(log((double) maxRank)/log(2.0)) + 2.0) + 1.0e-7);
        
        vector<double> results; results.resize(maxRank, 0.0);
        int nSamples = sampling.size();
        
        if (nSamples == 0) {  return results; }
        
        printf("M_PI = %f\n",M_PI);
        
#ifdef USE_GSL
        
        DiversityUtils dutils("lnabund");
        printf("nMax = %d \n",maxRank);
        for(int i = 0; i < sampling.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            printf("dmDash = %f dV = %f \n",sampling[i].alpha, sampling[i].beta);
            
            for (int j = 1; j <= maxRank; j++) {
                int nA = j;
                double dLog = 0.0, dP = 0.0;
                
                if(nA < 100){ //MAX_QUAD
                    dLog = logLikelihoodQuad2(nA, sampling[i].alpha, sampling[i].beta); //nA, dMDash, dV
                }
                else{
                    dLog = logLikelihoodRampal2(nA, sampling[i].alpha, sampling[i].beta); //nA, dMDash, dV
                }
                
                
                dP = exp(dLog);
                
                results[j - 1] += dP*sampling[i].ns;
                cout << j << '\t' << dP << '\t' << results[j-1] << endl;
            }
            
            
            
        }
        
        for (int i = 1; i<=maxRank; i++) {
            results[i-1] /= (double)nSamples;
            
            if (isnan(results[i-1]) || isinf(results[i-1])) { results[i-1] = 0.0; }
            
            cout <<i << '\t' << results[i-1] << endl;
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


