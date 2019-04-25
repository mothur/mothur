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

#define INIT_A            1.0
#define INIT_B            5.0
#define INIT_SIMPLEX_SIZE 0.1
#define INIT_S_SS         0.1
#define MIN_SIMPLEX_SIZE  1.0e-2
#define MAX_SIMPLEX_ITER  100000
#define MAX_LINE_LENGTH   1024
#define SAMPLE_FILE_SUFFIX ".sample"

#define DEF_SIGMA     0.1
#define DEF_SIGMA_S   100.0

#define DEF_ITER   250000
#define SLICE      10
#define PRECISION  1.0e-10

#ifdef USE_GSL

/***********************************************************************/
double chao(t_Data *ptData)
{
    double n1 = 0.0, n2 = 0.0;
    int **aanAbund = ptData->aanAbund;
    
    if(aanAbund[0][0] == 1 && aanAbund[1][0] == 2){
        n1 = (double) aanAbund[0][1]; n2 = (double) aanAbund[1][1];
        
        return ((double) ptData->nL) + 0.5*((n1*n1)/n2);
    }
    else{
        return -1.0;
    }
}
/***********************************************************************/
double nLogLikelihood(const gsl_vector * x, void * params)
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
int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params))
{
    const gsl_multimin_fminimizer_type *T =
    gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss;
    gsl_multimin_function minex_func;
    size_t iter = 0;
    int i = 0, status;
    double size;
    
    /* Initial vertex size vector */
    ss = gsl_vector_alloc (nP);
    
    /* Set all step sizes to default constant */
    gsl_vector_set_all(ss, INIT_SIMPLEX_SIZE);
    
    gsl_vector_set(ss,nP - 1,INIT_S_SS*gsl_vector_get(ptX,0));
    
    /* Initialize method and iterate */
    minex_func.f = f;
    minex_func.n = nP;
    minex_func.params = pvData;
    
    s = gsl_multimin_fminimizer_alloc (T, nP);
    gsl_multimin_fminimizer_set(s, &minex_func, ptX, ss);
    
    do{
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if(status)
            break;
        
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, MIN_SIMPLEX_SIZE);
        
        if(status == GSL_SUCCESS){
            for(i = 0; i < nP; i++){
                gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
            }
        }
    }
    while(status == GSL_CONTINUE && iter < MAX_SIMPLEX_ITER);
    
    for(i = 0; i < nP; i++){
        gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
    }
    
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    return status;
}
/***********************************************************************/
void outputResults(gsl_vector *ptX, t_Data *ptData)
{
    double dAlpha = 0.0, dBeta = 0.0, dS = 0.0, dL = 0.0;
    
    dAlpha = gsl_vector_get(ptX, 0);
    
    dBeta  = gsl_vector_get(ptX, 1);
    
    dS = gsl_vector_get(ptX, 2);
    
    dL = nLogLikelihood(ptX, ptData);
    
    MothurOut* m;  m = MothurOut::getInstance();
    m->mothurOut("\nML simplex: a = " + toString(dAlpha) +  " b = " + toString(dBeta) +  " S = " + toString(dS) +  " NLL = " + toString(dL) + "\n");
}
/***********************************************************************/
void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams)
{
    double dDeltaS =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaS);
    double dDeltaA =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaA);
    double dDeltaB =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaB);
    int    nSDash = 0;
    
    gsl_vector_set(ptXDash, 0, gsl_vector_get(ptX,0) + dDeltaA);
    gsl_vector_set(ptXDash, 1, gsl_vector_get(ptX,1) + dDeltaB);
    
    nSDash = nS + (int) floor(dDeltaS);
    if(nSDash < 1){
        nSDash = 1;
    }
    (*pnSDash) = nSDash;
}
/***********************************************************************/
double negLogLikelihood(double dAlpha, double dBeta, int nS, void * params)
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
void* metropolis (void * pvInitMetro)
{
    t_MetroInit *ptMetroInit  = (t_MetroInit *) pvInitMetro;
    gsl_vector  *ptX          = ptMetroInit->ptX;
    t_Data      *ptData       = ptMetroInit->ptData;
    t_Params    *ptParams     = ptMetroInit->ptParams;
    gsl_vector  *ptXDash      = gsl_vector_alloc(3); /*proposal*/
    char *szSampleFile = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));
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
    
    dNLL = negLogLikelihood(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), nS,(void*) ptData);
    
    string filename = ptParams->szOutFileStub + "_" + toString(ptMetroInit->nThread) + SAMPLE_FILE_SUFFIX;
    
    ofstream out; Utils util; util.openOutputFile(filename, out);
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        double dA = 0.0, dNLLDash = 0.0;
        
        getProposal(ptGSLRNG, ptXDash, ptX, &nSDash, nS, ptParams);
        
        dNLLDash = negLogLikelihood(gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), nSDash, (void*) ptData);
        
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

/***********************************************************************/
void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX)
{
    pthread_t thread1, thread2, thread3;
    int       iret1  , iret2  , iret3;
    gsl_vector *ptX1 = gsl_vector_alloc(3),
    *ptX2 = gsl_vector_alloc(3),
    *ptX3 = gsl_vector_alloc(3);
    t_MetroInit atMetroInit[3];
   
    
    MothurOut* m;  m = MothurOut::getInstance();
    m->mothurOut("\nMCMC iter = " + toString(ptParams->nIter) + " sigmaA = " + toString(ptParams->dSigmaA) +  " sigmaB = " + toString(ptParams->dSigmaB) +  " sigmaS = " + toString(ptParams->dSigmaS) + "\n");
    
    gsl_vector_memcpy(ptX1, ptX);
    
    gsl_vector_set(ptX2, 0, gsl_vector_get(ptX,0) + 2.0*ptParams->dSigmaA);
    gsl_vector_set(ptX2, 1, gsl_vector_get(ptX,1) + 2.0*ptParams->dSigmaB);
    gsl_vector_set(ptX2, 2, gsl_vector_get(ptX,2) + 2.0*ptParams->dSigmaS);
    
    gsl_vector_set(ptX3, 0, gsl_vector_get(ptX,0) - 2.0*ptParams->dSigmaA);
    gsl_vector_set(ptX3, 1, gsl_vector_get(ptX,1) - 2.0*ptParams->dSigmaB);
    if(gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaS > (double) ptData->nL){
        gsl_vector_set(ptX3, 2, gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaS);
    }
    else{
        gsl_vector_set(ptX3, 2, (double) ptData->nL);
    }
    atMetroInit[0].ptParams = ptParams;
    atMetroInit[0].ptData   = ptData;
    atMetroInit[0].ptX      = ptX1;
    atMetroInit[0].nThread  = 0;
    atMetroInit[0].lSeed    = ptParams->lSeed;
    atMetroInit[0].nAccepted = 0;
    
    //write thread 0
    m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX1, 0)) +  " b = " + toString(gsl_vector_get(ptX1, 1)) +  " S = " + toString(gsl_vector_get(ptX1, 2)) + "\n");

    
    atMetroInit[1].ptParams = ptParams;
    atMetroInit[1].ptData   = ptData;
    atMetroInit[1].ptX      = ptX2;
    atMetroInit[1].nThread  = 1;
    atMetroInit[1].lSeed    = ptParams->lSeed + 1;
    atMetroInit[1].nAccepted = 0;
    
    //write thread 1
    m->mothurOut(toString(atMetroInit[1].nThread) + ": a = " + toString(gsl_vector_get(ptX2, 0)) +  " b = " + toString(gsl_vector_get(ptX2, 1)) +  " S = " + toString(gsl_vector_get(ptX2, 2)) + "\n");
    
    
    atMetroInit[2].ptParams = ptParams;
    atMetroInit[2].ptData   = ptData;
    atMetroInit[2].ptX      = ptX3;
    atMetroInit[2].nThread  = 2;
    atMetroInit[2].lSeed    = ptParams->lSeed + 2;
    atMetroInit[2].nAccepted = 0;
    
    //write thread 2
    m->mothurOut(toString(atMetroInit[2].nThread) + ": a = " + toString(gsl_vector_get(ptX3, 0)) +  " b = " + toString(gsl_vector_get(ptX3, 1)) +  " S = " + toString(gsl_vector_get(ptX3, 2)) + "\n");

    
    iret1 = pthread_create(&thread1, NULL, metropolis, (void*) &atMetroInit[0]);
    iret2 = pthread_create(&thread2, NULL, metropolis, (void*) &atMetroInit[1]);
    iret3 = pthread_create(&thread3, NULL, metropolis, (void*) &atMetroInit[2]);
    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);
    pthread_join(thread3, NULL);
    
    
    m->mothurOut(toString(atMetroInit[0].nThread) +": accept. ratio " + toString(atMetroInit[0].nAccepted) + "/" + toString(ptParams->nIter) +  " = " + toString(((double) atMetroInit[0].nAccepted)/((double) ptParams->nIter)) +  "\n");
    m->mothurOut(toString(atMetroInit[1].nThread) +": accept. ratio " + toString(atMetroInit[1].nAccepted) + "/" + toString(ptParams->nIter) +  " = " + toString(((double) atMetroInit[1].nAccepted)/((double) ptParams->nIter)) +  "\n");
    m->mothurOut(toString(atMetroInit[2].nThread) +": accept. ratio " + toString(atMetroInit[2].nAccepted) + "/" + toString(ptParams->nIter) +  " = " + toString(((double) atMetroInit[2].nAccepted)/((double) ptParams->nIter)) +  "\n");
    
    gsl_vector_free(ptX1); gsl_vector_free(ptX2); gsl_vector_free(ptX3);
}

#endif
/***********************************************************************/
void MetroIG::loadAbundance(t_Data *ptData, SAbundVector* rank)
{
    int nNA = 0; int  nL = 0, nJ = 0;
    int maxRank = rank->getMaxRank();
    
    for(int i = 1; i <= maxRank; i++){ if (rank->get(i) != 0) { nNA++; } }
    
    int **aanAbund = NULL;
    aanAbund = (int **) malloc(nNA*sizeof(int*));
    
    int count = 0;
    for(int i = 1; i <= maxRank; i++){
        
        if (rank->get(i) != 0) {
            aanAbund[count] = (int *) malloc(sizeof(int)*2);
            
            int nA = i;
            int nC = rank->get(i);
            
            nL += nC;
            nJ += nC*nA;
            
            aanAbund[count][0]  = nA;
            aanAbund[count][1]  = nC;
            
            count++;
        }
        
    }
    
    ptData->nJ          = nJ;
    ptData->nL          = nL;
    ptData->aanAbund    = aanAbund;
    ptData->nNA         = nNA;
}


/***********************************************************************/
void MetroIG::freeAbundance(t_Data *ptData)
{
    for(int i = 0; i < ptData->nNA; i++){
        free(ptData->aanAbund[i]);
    }
    free(ptData->aanAbund);
}

/***********************************************************************/
vector<string> MetroIG::getValues(SAbundVector* rank){
    try {
    
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaA = sigmaA; tParams.dSigmaB = sigmaB; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub; tParams.lSeed = m->getRandomSeed();
        t_Data   tData;
#ifdef USE_GSL
        loadAbundance(&tData, rank);

        gsl_vector* ptX = gsl_vector_alloc(3); /*parameter estimates*/
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl

        gsl_rng_env_setup();
        
        gsl_set_error_handler_off();

        /*set initial estimates for parameters*/
        gsl_vector_set(ptX, 0, INIT_A);
        gsl_vector_set(ptX, 1, INIT_B);
        gsl_vector_set(ptX, 2, numOTUs*2);
        
        double chaoResult = chao(&tData);
        m->mothurOut("\nD = " + toString(numOTUs) + " L = " + toString(sampled) +  " Chao = " + toString(chaoResult) +  "\n");
        
        minimiseSimplex(ptX, 3, (void*) &tData, &nLogLikelihood);
        
        outputResults(ptX, &tData);
        
        if(tParams.nIter > 0){
           mcmc(&tParams, &tData, ptX);
        }
        
        /*free up allocated memory*/
        gsl_vector_free(ptX);
      
        freeAbundance(&tData);
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
