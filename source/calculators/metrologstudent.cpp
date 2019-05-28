//
//  metrologstudent.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/2/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metrologstudent.hpp"

/***********************************************************************/
MetroLogStudent::MetroLogStudent(double sigm, double sigv, double sign, double sigS, int n, string st) : sigmaM(sigm), sigmaV(sigv), sigmaN(sign), sigmaS(sigS), nIters(n), outFileStub(st), DiversityCalculator(false) {}
/***********************************************************************/


#ifdef USE_GSL
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
    
    DiversityUtils dutils("metrols");
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < 100){ //MAX_QUAD
            dLogP = dutils.logLikelihoodQuad(nA, dMDash, dV, dNu);
        }
        else{
            dLogP = dutils.logLikelihoodRampal(nA, dMDash, dV, dNu);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = dutils.logLikelihoodQuad(0, dMDash, dV, dNu);
    
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
    DiversityUtils dutils("metrols");
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < 100){ //MAX_QUAD
            dLogP = dutils.logLikelihoodQuad(nA, dMDash, dV, dNu);
        }
        else{
            dLogP = dutils.logLikelihoodRampal(nA, dMDash, dV, dNu);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = dutils.logLikelihoodQuad(0, dMDash, dV, dNu);
    
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
        }else if(nIter % 10000 == 0){
            cout << nIter << endl;
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
#endif
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
