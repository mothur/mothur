//
//  metrologstudent.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/2/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metrologstudent.hpp"

/***********************************************************************/
double f1_1(double x, void *pvParams)
{
    t_LSParams *ptLSParams = (t_LSParams *) pvParams;
    double dMDash = ptLSParams->dMDash, dV = ptLSParams->dV, dNu = ptLSParams->dNu;
    int n = ptLSParams->n;
    double t = ((x - dMDash)*(x - dMDash))/dV;
    double dExp  = x*((double) n) - exp(x);
    double dF  = pow(1.0 + t/dNu, -0.5*(dNu + 1.0));
    
    return exp(dExp)*dF;
}
/***********************************************************************/
double f1Log(double x, void *pvParams)
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
double logStirlingsGamma(double dZ)
{
    return 0.5*log(2.0*M_PI) + (dZ - 0.5)*log(dZ) - dZ;
}
/***********************************************************************/
double logLikelihoodQuad(int n, double dMDash, double dV, double dNu)
{
    gsl_integration_workspace *ptGSLWS =
    gsl_integration_workspace_alloc(1000);
    double dLogFac1   = 0.0, dLogFacN  = 0.0;
    double dN = (double) n, dResult = 0.0, dError = 0.0, dPrecision = 0.0;
    gsl_function tGSLF;
    t_LSParams tLSParams;
    double dA = 0.0, dB = 0.0;
    
    tLSParams.n = n; tLSParams.dMDash = dMDash; tLSParams.dV = dV; tLSParams.dNu = dNu;
    
    tGSLF.function = &f1_1;
    tGSLF.params   = (void *) &tLSParams;
    
    if(dNu < 100){ //MAX_MU_GAMMA
        dLogFac1 = gsl_sf_lngamma(0.5*(dNu + 1.0)) - gsl_sf_lngamma(0.5*dNu) - 0.5*log(M_PI*dNu);
    }
    else{
        dLogFac1 = 0.5*dNu*(log(0.5*(dNu + 1.0)) - log(0.5*dNu)) -0.5*log(2.0*M_PI) - 0.5;
    }
    
    if(n < 50){
        dLogFacN = gsl_sf_fact(n);
        dLogFacN = log(dLogFacN);
    }
    else if(n < 100){
        dLogFacN = gsl_sf_lngamma(dN + 1.0);
    }
    else{
        dLogFacN = logStirlingsGamma(dN + 1.0);
    }
    
    dA = -100.0; dB = 100.0;
    
    if(n < 10){
        dPrecision = HI_PRECISION;
    }
    else{
        dPrecision = LO_PRECISION;
    }
    
    gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
    
    //printf("%f %f\n", dResult, dError);
    
    gsl_integration_workspace_free(ptGSLWS);
    
    return log(dResult) - dLogFacN + dLogFac1 - 0.5*log(dV);
}
/***********************************************************************/
double logLikelihoodLNQuad(int n, double dMDash, double dV)
{
    gsl_integration_workspace *ptGSLWS =
    gsl_integration_workspace_alloc(1000);
    double dLogFac1   = 0.0, dLogFacN  = 0.0;
    double dResult = 0.0, dError = 0.0, dPrecision = 0.0;
    gsl_function tGSLF;
    t_LNParams tLNParams;
    double dEst = dMDash + ((double)n)*dV, dA = 0.0, dB = 0.0;
    
    tLNParams.n = n; tLNParams.dMDash = dMDash; tLNParams.dV = dV;
    
    tGSLF.function = &f1Log;
    tGSLF.params   = (void *) &tLNParams;
    
    dLogFac1 = log(2.0*M_PI*dV);
    
    if(n < 50){
        dLogFacN = gsl_sf_fact(n);
        dLogFacN = log(dLogFacN);
    }
    else{
        dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
    }
    
    DiversityUtils dutils("metrols");
    
    if(dEst > dV){
        double dMax = 0.0;
        double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
        double dVar   = 0.0;
        
        if(fabs(dUpper) > 1.0e-7){
            dutils.solveF(0.0, dUpper, (void *) &tLNParams, 1.0e-5, &dMax);
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
            dutils.solveF(dLower, dUpper,  (void *) &tLNParams, 1.0e-5, &dMax);
        }
        else{
            dMax = 0.5*(dLower + dUpper);
        }
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
double logLikelihoodLNRampal(int n, double dMDash, double dV)
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
double logLikelihoodRampal(int n, double dMDash, double dV, double dNu)
{
    double dGamma = 0.5*(dNu + 1.0), dN = (double) n, dRN = 1.0/dN, dRSV = 1.0/(sqrt(dV)*sqrt(dNu));
    double dZ = (log(dN) - dMDash)*dRSV;
    double dDZDX = dRN*dRSV, dDZDX2 = -dRN*dRN*dRSV;
    double dF = (1.0 + dZ*dZ);
    double dA = 0.0, dB = 0.0, dTemp = 0.0;
    double dLogFac1 = 0.0;
    
    if(dNu < 100){ //MAX_MU_GAMMA
        dLogFac1 = gsl_sf_lngamma(0.5*(dNu + 1.0)) - gsl_sf_lngamma(0.5*dNu) - 0.5*log(M_PI*dNu);
    }
    else{
        dLogFac1 = 0.5*dNu*(log(0.5*(dNu + 1.0)) - log(0.5*dNu)) -0.5*log(2.0*M_PI) - 0.5;
    }
    
    dA = 4.0*dZ*dZ*dDZDX*dDZDX*dGamma*(dGamma + 1.0);
    dA /= dF*dF;
    
    dB = -2.0*dGamma*(dDZDX*dDZDX + dZ*dDZDX2);
    dB /= dF;
    
    dTemp = dRN + dA + dB;
    
    return -dGamma*log(dF) + log(dTemp) + dLogFac1 - 0.5*log(dV);
}

/***********************************************************************/
double nLogLikelihood2(const gsl_vector * x, void * params)
{
    double dMDash  = gsl_vector_get(x,0), dV = gsl_vector_get(x,1);
    double dNu  = gsl_vector_get(x,2);
    int    nS = (int) floor(gsl_vector_get(x, 3));
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogNot0 = 0.0, dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dV <= 0.0 || dNu < 1.0){
        return PENALTY;
    }
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < 100){ //MAX_QUAD
            dLogP = logLikelihoodQuad(nA, dMDash, dV, dNu);
        }
        else{
            dLogP = logLikelihoodRampal(nA, dMDash, dV, dNu);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = logLikelihoodQuad(0, dMDash, dV, dNu);
    
    dLog1 = (nS - ptData->nL)*dLog0;
    
    dLog2 = - gsl_sf_lnfact(nS - ptData->nL);
    
    dLog3 = gsl_sf_lnfact(nS);
    
    dLogL += dLog1 + dLog2 + dLog3;
    
    /*return*/
    return -dLogL;
}
/***********************************************************************/
double negLogLikelihood(double dMDash, double dV, double dNu, int nS, void * params)
{
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dV <= 0.0 || dNu < 1.0){
        return PENALTY;
    }
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < 100){ //MAX_QUAD
            dLogP = logLikelihoodQuad(nA, dMDash, dV, dNu);
        }
        else{
            dLogP = logLikelihoodRampal(nA, dMDash, dV, dNu);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = logLikelihoodQuad(0, dMDash, dV, dNu);
    
    dLog1 = (nS - ptData->nL)*dLog0;
    
    dLog2 = - gsl_sf_lnfact(nS - ptData->nL);
    
    dLog3 = gsl_sf_lnfact(nS);
    
    dLogL += dLog1 + dLog2 + dLog3;
    
    /*return*/
    return -dLogL;
}
/***********************************************************************/
void* metropolis2 (void * pvInitMetro)
{
    t_MetroInit *ptMetroInit  = (t_MetroInit *) pvInitMetro;
    gsl_vector  *ptX          = ptMetroInit->ptX;
    t_Data      *ptData       = ptMetroInit->ptData;
    t_Params    *ptParams     = ptMetroInit->ptParams;
    gsl_vector  *ptXDash      = gsl_vector_alloc(4); /*proposal*/
    char *szSampleFile = (char *) malloc(1024*sizeof(char));
    const gsl_rng_type *T;
    gsl_rng            *ptGSLRNG;
    int nS = 0, nSDash = 0, nIter = 0;
    double dRand = 0.0, dNLL = 0.0;
    void   *pvRet = NULL;
    
    /*set up random number generator*/
    T        = gsl_rng_default;
    ptGSLRNG = gsl_rng_alloc (T);
    
    nS = (int) floor(gsl_vector_get(ptX,3));
    
    dNLL = negLogLikelihood(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), gsl_vector_get(ptX,2), nS,(void*) ptData);
    
    string filename = ptParams->szOutFileStub + "_" + toString(ptMetroInit->nThread) + ".sample";
    
    ofstream out; Utils util; util.openOutputFile(filename, out);
    out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    DiversityUtils dutils("metrols");
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        double dA = 0.0, dNLLDash = 0.0;
        
        dutils.getProposal(ptGSLRNG, ptXDash, ptX, &nSDash, nS, ptParams);
        
        dNLLDash = negLogLikelihood(gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), gsl_vector_get(ptXDash,2), nSDash, (void*) ptData);
        
        dA = exp(dNLL - dNLLDash);
        
        if(dA > 1.0){ dA = 1.0; }
        
        dRand = gsl_rng_uniform(ptGSLRNG);
        
        if(dRand < dA){
            gsl_vector_memcpy(ptX, ptXDash);
            nS = nSDash;
            dNLL = dNLLDash;
            ptMetroInit->nAccepted++;
        }
        
        if(nIter % 10 == 0){
            out << nIter << "," << gsl_vector_get(ptX, 0) << "," << gsl_vector_get(ptX, 1) << "," << gsl_vector_get(ptX, 2) << "," << nS << "," << dNLL << endl;
        }
        
        nIter++;
    }
    
    out.close();
    
    /*free up allocated memory*/
    gsl_vector_free(ptXDash);
    free(szSampleFile);
    gsl_rng_free(ptGSLRNG);
    
    return pvRet;
}
/***********************************************************************/
vector<string> MetroLogStudent::getValues(SAbundVector* rank){
    try {
        
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaX = sigmaM; tParams.dSigmaY = sigmaV; tParams.dSigmaN = sigmaN; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub; tParams.lSeed = m->getRandomSeed();
        t_Data   tData;
        
#ifdef USE_GSL
        
        DiversityUtils dutils("metrols");
        
        dutils.loadAbundance(&tData, rank);
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl
        
        gsl_vector* ptX = gsl_vector_alloc(4); /*parameter estimates*/
        
        gsl_rng_env_setup();
        
        gsl_set_error_handler_off();
        
        /*set initial estimates for parameters*/
        gsl_vector_set(ptX, 0, -10.0); //INIT_M
        gsl_vector_set(ptX, 1, 20.0); //INIT_V
        gsl_vector_set(ptX, 2, 20.0); //INIT_N
        gsl_vector_set(ptX, 3, numOTUs*2);
        
        double chaoResult = dutils.chao(&tData);
        m->mothurOut("\nMetroLogStudent - D = " + toString(numOTUs) + " L = " + toString(sampled) +  " Chao = " + toString(chaoResult) +  "\n");
        
        dutils.minimiseSimplex(ptX, 4, (void*) &tData, &nLogLikelihood2, 0.1, 1.0e-3, 100000);
        dutils.outputResults(ptX, &tData, &nLogLikelihood2);
        
        if(tParams.nIter > 0){ dutils.mcmc(&tParams, &tData, ptX, &metropolis2); }
        
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
        m->errorOut(e, "MetroLogStudent", "getValues");
        exit(1);
    }
}
/***********************************************************************/
