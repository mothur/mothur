//
//  metroig.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metroig.hpp"

/*constants for simplex minimisation*/
#define PENALTY           1.0e20
#define SLICE      10


#ifdef USE_GSL


/***********************************************************************/
double nLogLikelihood0(const gsl_vector * x, void * params)
{
    double dAlpha  = gsl_vector_get(x,0), dBeta = gsl_vector_get(x,1);
    int    nS = (int) floor(gsl_vector_get(x, 2));
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dAlpha <= 0.0 || dBeta <= 0.0){
        return PENALTY;
    }
    
    DiversityUtils dutils;
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        dLogP = dutils.logLikelihood(nA, dAlpha, dBeta);
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = dutils.logLikelihood(0, dAlpha, dBeta);
    
    dLog1 = (nS - ptData->nL)*dLog0;
    
    dLog2 = - gsl_sf_lnfact(nS - ptData->nL);
    
    dLog3 = gsl_sf_lnfact(nS);
    
    dLogL += dLog1 + dLog2 + dLog3;
    
    /*return*/
    return -dLogL;
}
/***********************************************************************/
double negLogLikelihood0(double dAlpha, double dBeta, int nS, void * params)
{
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dAlpha <= 0.0 || dBeta <= 0.0){
        return PENALTY;
    }
    
    DiversityUtils dutils;
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        dLogP = dutils.logLikelihood(nA, dAlpha, dBeta);
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = dutils.logLikelihood(0, dAlpha, dBeta);
    
    dLog1 = (nS - ptData->nL)*dLog0;
    
    dLog2 = - gsl_sf_lnfact(nS - ptData->nL);
    
    dLog3 = gsl_sf_lnfact(nS);
    
    dLogL += dLog1 + dLog2 + dLog3;
    
    /*return*/
    return -dLogL;
}

/***********************************************************************/
void* metropolis0 (void * pvInitMetro)
{
    t_MetroInit *ptMetroInit  = (t_MetroInit *) pvInitMetro;
    gsl_vector  *ptX          = ptMetroInit->ptX;
    t_Data      *ptData       = ptMetroInit->ptData;
    t_Params    *ptParams     = ptMetroInit->ptParams;
    gsl_vector  *ptXDash      = gsl_vector_alloc(3); /*proposal*/
    char *szSampleFile = (char *) malloc(1024*sizeof(char));
    const gsl_rng_type *T;
    gsl_rng            *ptGSLRNG;
    //FILE    *sfp = NULL;
    int nS = 0, nSDash = 0, nIter = 0;
    double dRand = 0.0, dNLL = 0.0;
    void   *pvRet = NULL;
    
    /*set up random number generator*/
    T        = gsl_rng_default;
    ptGSLRNG = gsl_rng_alloc (T);
    
    nS = (int) floor(gsl_vector_get(ptX,2));
    
    dNLL = negLogLikelihood0(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), nS,(void*) ptData);
    
    string filename = ptParams->szOutFileStub + "_" + toString(ptMetroInit->nThread) + ".sample";
    
    ofstream out; Utils util; util.openOutputFile(filename, out);
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    DiversityUtils dutils;
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        double dA = 0.0, dNLLDash = 0.0;
        
        dutils.getProposal(ptGSLRNG, ptXDash, ptX, &nSDash, nS, ptParams);
        
        dNLLDash = negLogLikelihood0(gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), nSDash, (void*) ptData);
        
        dA = exp(dNLL - dNLLDash);
        if(dA > 1.0){
            dA = 1.0;
        }
        
        dRand = gsl_rng_uniform(ptGSLRNG);
        
        if(dRand < dA){
            gsl_vector_memcpy(ptX, ptXDash);
            nS = nSDash;
            dNLL = dNLLDash;
            ptMetroInit->nAccepted++;
        }
        
        if(nIter % SLICE == 0){
            out << nIter << "," << gsl_vector_get(ptX, 0) << "," << gsl_vector_get(ptX, 1) << "," << nS << "," << dNLL << endl;
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
vector<string> MetroIG::getValues(SAbundVector* rank){
    try {
        
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaX = sigmaA; tParams.dSigmaY = sigmaB; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub; tParams.lSeed = m->getRandomSeed();
        t_Data   tData;
#ifdef USE_GSL
        
        DiversityUtils dutils;
        
        dutils.loadAbundance(&tData, rank);
        
        gsl_vector* ptX = gsl_vector_alloc(3); /*parameter estimates*/
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl
        
        gsl_rng_env_setup();
        
        gsl_set_error_handler_off();
        
        /*set initial estimates for parameters*/
        gsl_vector_set(ptX, 0, 1.0);
        gsl_vector_set(ptX, 1, 5.0);
        gsl_vector_set(ptX, 2, numOTUs*2);
        
        double chaoResult = dutils.chao(&tData);
        m->mothurOut("\nMetroIG - D = " + toString(numOTUs) + " L = " + toString(sampled) +  " Chao = " + toString(chaoResult) +  "\n");
        
        dutils.minimiseSimplex(ptX, 3, (void*) &tData, &nLogLikelihood0, 0.1, "metroig");
        
        dutils.outputResults(ptX, &tData, &nLogLikelihood0, "metroig");
        
        if(tParams.nIter > 0){ dutils.mcmc(&tParams, &tData, ptX, &metropolis0); }
        
        /*free up allocated memory*/
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
        m->errorOut(e, "MetroIG", "getValues");
        exit(1);
    }
}
/***********************************************************************/

