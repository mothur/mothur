//
//  metrolognormal.c
//  Mothur
//
//  Created by Sarah Westcott on 4/25/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metrolognormal.hpp"

#define HI_PRECISION 1.0e-12
#define LO_PRECISION 1.0e-7

/*constants for calculated compound Poisson lognormal*/
#define V_MULT          25.0
#define MAX_QUAD        100
#define MAX_QUAD_DERIV  100

/***********************************************************************/
double f1(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV;
    int n = ptLNParams->n;
    double dTemp = (x - dMDash);
    double dExp  = x*((double) n) - exp(x) - 0.5*((dTemp*dTemp)/dV);
    double dRet  = exp(dExp);
    
    return dRet;
}
/***********************************************************************/
double derivExponent(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
    double dTemp = (x - dMDash)/dV, dRet = 0.0;
    
    dRet = ((double) n) - exp(x) - dTemp;
    
    return dRet;
}
/***********************************************************************/
double logLikelihoodRampal(int n, double dMDash, double dV)
{
    double dN = (double) n;
    double dLogLik = 0.0, dTemp = gsl_pow_int(log(dN) - dMDash,2), dTemp3 = gsl_pow_int(log(dN) - dMDash,3);
    
    dLogLik = -0.5*log(2.0*M_PI*dV) - log(dN) - (dTemp/(2.0*dV));
    
    dLogLik += log(1.0 + 1.0/(2.0*dN*dV)*(dTemp/dV + log(dN) - dMDash - 1.0)
                   + 1.0/(6.0*dN*dN*dV*dV*dV)*(3.0*dV*dV - (3.0*dV - 2.0*dV*dV)*(dMDash - log(dN))
                                               - 3.0*dV*dTemp + dTemp3));
    
    return dLogLik;
}
/***********************************************************************/
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
    
    //printf("%f %f %d %f %f\n",ptLNParams->dMDash, ptLNParams->dV, ptLNParams->n, x_lo, x_hi);
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
    
    tGSLF.function = &f1;
    tGSLF.params   = (void *) &tLNParams;
    
    dLogFac1 = log(2.0*M_PI*dV);
    
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
        
        if(fabs(dUpper) > 1.0e-7){
            solveF(0.0, dUpper, derivExponent, (void *) &tLNParams, 1.0e-5, &dMax);
        }
        
        dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
        
        dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
    }
    else{
        double dMax = 0.0;
        double dLower = dEst - dV;
        double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
        double dVar   = 0.0;
        
        if(fabs(dUpper - dLower) > 1.0e-7){
            solveF(dLower, dUpper, derivExponent, (void *) &tLNParams, 1.0e-5, &dMax);
        }
        else{
            dMax = 0.5*(dLower + dUpper);
        }
        dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
        
        dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
    }
    
    if(n < 10)  {  dPrecision = HI_PRECISION;   }
    else        {  dPrecision = LO_PRECISION;   }
    
    gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
    
    gsl_integration_workspace_free(ptGSLWS);
    
    return log(dResult) - dLogFacN -0.5*dLogFac1;
}
/***********************************************************************/
double nLogLikelihood1(const gsl_vector * x, void * params)
{
    double dMDash  = gsl_vector_get(x,0), dV = gsl_vector_get(x,1);
    int    nS = (int) floor(gsl_vector_get(x, 2));
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogL   = 0.0;
    
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < MAX_QUAD){
            dLogP = logLikelihoodQuad(nA, dMDash, dV);
        }
        else{
            dLogP = logLikelihoodRampal(nA, dMDash, dV);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
    }
    
    dLogL += (nS - ptData->nL)*logLikelihoodQuad(0, dMDash, dV);
    
    dLogL -= gsl_sf_lnfact(nS - ptData->nL);
    
    dLogL += gsl_sf_lnfact(nS);
   
    /*return*/
    return -dLogL;
}
/***********************************************************************/
double negLogLikelihood1(double dMDash, double dV, int nS, void * params)
{
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLog0 = 0.0, dLogL   = 0.0;
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < MAX_QUAD){
            dLogP = logLikelihoodQuad(nA, dMDash, dV);
        }
        else{
            dLogP = logLikelihoodRampal(nA, dMDash, dV);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
    }
    
    dLog0     = logLikelihoodQuad(0, dMDash, dV);
    
    if(nS > ptData->nL){
        dLogL += (nS - ptData->nL)*dLog0;
    }
    
    dLogL -= gsl_sf_lnfact(nS - ptData->nL);
    
    dLogL += gsl_sf_lnfact(nS);
    
    /*return*/
    return -dLogL;
}
/***********************************************************************/
void getProposal1(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX,
                 int* pnSDash, int nS, t_MetroInit *ptMetroInit)
{
    t_Params *ptParams = ptMetroInit->ptParams;
    double    dDelta1 = gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaX);
    double    dDelta2 = gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaY);
    double    dDeltaS = gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaS);
    int       nSDash = 0;
    
    gsl_vector_set(ptXDash, 0, gsl_vector_get(ptX,0) + dDelta1);
    gsl_vector_set(ptXDash, 1, gsl_vector_get(ptX,1) + dDelta2);
    
    nSDash = nS + (int) floor(dDeltaS);
    if(nSDash < 1){
        nSDash = 1;
    }
    (*pnSDash) = nSDash;
}
/***********************************************************************/
void* metropolis1 (void * pvInitMetro)
{
    t_MetroInit *ptMetroInit  = (t_MetroInit *) pvInitMetro;
    gsl_vector  *ptX          = ptMetroInit->ptX;
    t_Data      *ptData       = ptMetroInit->ptData;
    t_Params    *ptParams     = ptMetroInit->ptParams;
    gsl_vector  *ptXDash      = gsl_vector_alloc(3); /*proposal*/
    
    const gsl_rng_type *T;
    gsl_rng            *ptGSLRNG;
    
    int nS = 0, nSDash = 0,nIter = 0;
    double dRand = 0.0, dNLL = 0.0;
    void   *pvRet = NULL;
    double dM = 0.0, dV = 0.0;
    double dMDash = 0.0, dVDash = 0.0;
    double dXDash = 0.0, dX = 0.0;
    
    /*set up random number generator*/
    T        = gsl_rng_default;
    ptGSLRNG = gsl_rng_alloc (T);
    
    nS = (int) floor(gsl_vector_get(ptX,2));
    
    dNLL = negLogLikelihood1(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), nS,(void*) ptData);
    
    dM = gsl_vector_get(ptX,0); dV = gsl_vector_get(ptX,1);
    gsl_vector_set(ptX,0,dM + 0.5*dV);
        
    string filename = ptParams->szOutFileStub + "_" + toString(ptMetroInit->nThread) + ".sample";
    
    ofstream out; Utils util; util.openOutputFile(filename, out);
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        double dA = 0.0, dNLLDash = 0.0;
        
        getProposal1(ptGSLRNG, ptXDash, ptX, &nSDash, nS,ptMetroInit);
        
        dXDash =  gsl_vector_get(ptXDash,0); dVDash = gsl_vector_get(ptXDash,1);
        dMDash =  dXDash - 0.5*dVDash;
        dNLLDash = negLogLikelihood1(dMDash, dVDash, nSDash, (void*) ptData);
        
        dA = exp(dNLL - dNLLDash);
        if(dA > 1.0){
            dA = 1.0;
        }
        
        dRand = gsl_rng_uniform(ptGSLRNG);
        
        if(dRand < dA){
            ptMetroInit->nAccepted++;
            gsl_vector_memcpy(ptX, ptXDash);
            nS = nSDash;
            dNLL = dNLLDash;
        }
        
        if(nIter % 10 == 0){
            dX =  gsl_vector_get(ptX,0); dV = gsl_vector_get(ptX,1);
            dM =  dX - 0.5*dV;
            
            out << nIter << "," << dM << "," << dV << "," << nS << "," << dNLL << endl;
            
        }
        
        nIter++;
    }
    out.close();
    
    /*free up allocated memory*/
    gsl_vector_free(ptXDash);
    gsl_rng_free(ptGSLRNG);
    
    return pvRet;
}
/***********************************************************************/
void MetroLogNormal::outputResults(gsl_vector *ptX, t_Data *ptData)
{
    double dMDash = 0.0, dV = 0.0, dS = 0.0, dL = 0.0;
    
    dMDash = gsl_vector_get(ptX, 0);
    
    dV     = gsl_vector_get(ptX, 1);
    
    dS     = gsl_vector_get(ptX, 2);
    
    dL = nLogLikelihood1(ptX, ptData);
    
    m->mothurOut("\nMetroLogNormal - ML simplex: M = " + toString(dMDash) +  " V = " + toString(dV) +  " S = " + toString(dS) +  " NLL = " + toString(dL) + "\n");
}
/***********************************************************************/
vector<string> MetroLogNormal::getValues(SAbundVector* rank){
    try {
        
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaX = sigmaX; tParams.dSigmaY = sigmaY; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub; tParams.lSeed = m->getRandomSeed();
        t_Data   tData;
        
#ifdef USE_GSL
        
        DiversityUtils dutils;
        
        dutils.loadAbundance(&tData, rank);

        gsl_vector* ptX = gsl_vector_alloc(3);
        
        gsl_rng_env_setup();

        gsl_set_error_handler_off();

        dutils.loadAbundance(&tData, rank);
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl

        gsl_vector_set(ptX, 0, INIT_M_DASH);
        gsl_vector_set(ptX, 1, INIT_V);
        gsl_vector_set(ptX, 2, numOTUs*2);

        double chaoResult = dutils.chao(&tData);
        m->mothurOut("\nMetroLogNormal - D = " + toString(numOTUs) + " L = " + toString(sampled) +  " Chao = " + toString(chaoResult) +  "\n");

        dutils.minimiseSimplex(ptX, 3, (void*) &tData, &nLogLikelihood1, 1.0, "metroln");

        outputResults(ptX, &tData);

        if(tParams.nIter > 0){ dutils.mcmc(&tParams, &tData, ptX, &metropolis1); }
       
        gsl_vector_free(ptX);

        dutils.freeAbundance(&tData);
#endif
        
        vector<string> outputs;
        outputs.push_back(outFileStub + "_0.sample");
        outputs.push_back(outFileStub + "_1.sample");
        outputs.push_back(outFileStub + "_2.sample");
        
        return outputs;
    }
    catch(exception& e) {
        m->errorOut(e, "MetroLogNormal", "getValues");
        exit(1);
    }
}
/***********************************************************************/


