//
//  metroig.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metroig.hpp"

/*User includes*/
#include "FileUtils.h"
#include "MatrixUtils.h"

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
#define TRUE  1
#define FALSE 0

//static int  nLines   = 11;

static int  verbose  = FALSE;

#ifdef USE_GSL
/***********************************************************************/
double fX(double x, double dA, double dB, double dNDash)
{
    double dTemp1 = (dA*(x - dB)*(x - dB))/x;
    
    return log(x) - (1.0/dNDash)*(x + dTemp1);
}
/***********************************************************************/
double f2X(double x, double dA, double dB, double dNDash)
{
    double dRet = 0.0, dTemp = 2.0*dA*dB*dB;
    
    dRet = (1.0/(x*x))*(1.0 + (1.0/dNDash)*(dTemp/x));
    
    return -dRet;
}
/***********************************************************************/
double sd(int n, double dAlpha, double dBeta)
{
    double dGamma = -0.5;
    double dA = 0.5*(-1.0 + sqrt(1.0 + (dAlpha*dAlpha)/(dBeta*dBeta)));
    double dN = (double) n, dNDash = dN + dGamma - 1.0, dRN = 1.0/dN;
    double dTemp1 = (0.5*dN)/(1.0 + dA), dTemp2 = 4.0*dRN*dRN*(1.0 + dA)*dA*dBeta*dBeta;
    double dXStar = dTemp1*(1.0 + sqrt(1.0 + dTemp2));
    double dFX = fX(dXStar, dA, dBeta, dNDash);
    double d2FX = -dNDash*f2X(dXStar, dA, dBeta, dNDash);
    double dLogK = 0.0, dGamma1 = dGamma;
    
    if(dGamma1 < 0.0){
        dGamma1 *= -1.0;
    }
    
    dLogK = gsl_sf_bessel_lnKnu(dGamma1,2.0*dA*dBeta);
    
    return -2.0*dA*dBeta -log(2.0) -dLogK -dGamma*log(dBeta) + dNDash*dFX + 0.5*log(2.0*M_PI) - 0.5*log(d2FX);
}
/***********************************************************************/
int bessel(double* pdResult, int n, double dAlpha, double dBeta)
{
    double dGamma  = -0.5;
    double dResult = 0.0;
    double dOmega = 0.0, dGamma2 = 0.0;
    double dLogK1 = 0.0, dLogK2 = 0.0;
    double dN = (double) n, dNu = dGamma + dN;
    double dTemp1 = 0.0;
    
    if(dNu < 0.0){
        dNu = -dNu;
    }
    
    if(dGamma < 0.0){
        dGamma2 = -dGamma;
    }
    else{
        dGamma2 = dGamma;
    }
    
    dOmega = sqrt(dBeta*dBeta + dAlpha*dAlpha) - dBeta;
    
    dLogK2 = gsl_sf_bessel_lnKnu(dNu, dAlpha);
    
    if(!gsl_finite(dLogK2)){
        if(dAlpha < 0.1*sqrt(dNu + 1.0)){
            //printf("l ");
            dLogK2 = gsl_sf_lngamma(dNu) + (dNu - 1.0)*log(2.0) - dNu*log(dAlpha);
        }
        else{
            //printf("sd ");
            (*pdResult) = dResult;
            return FALSE;
        }
    }
    
    dLogK1 = dGamma*log(dOmega/dAlpha) -gsl_sf_bessel_lnKnu(dGamma2,dOmega);
    
    dTemp1 = log((dBeta*dOmega)/dAlpha);
    
    dResult = dN*dTemp1 + dLogK2 + dLogK1;
    (*pdResult) = dResult;
    return TRUE;
}

/***********************************************************************/
double logLikelihood(int n, double dAlpha, double dBeta)
{
    double dLogFacN = 0.0;
    int status      = 0;
    double dRet     = 0.0;
    
    if(n < 50){
        dLogFacN = gsl_sf_fact(n);
        dLogFacN = log(dLogFacN);
    }
    else{
        dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
    }
    
    status = bessel(&dRet,n, dAlpha,dBeta);
    if(status == FALSE){
        dRet = sd(n, dAlpha,dBeta);
    }
    
    return dRet - dLogFacN;
}
/***********************************************************************/
double nLogLikelihood(const gsl_vector * x, void * params)
{
    double dAlpha  = gsl_vector_get(x,0), dBeta = gsl_vector_get(x,1);
    int    nS = (int) floor(gsl_vector_get(x, 2));
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogNot0 = 0.0, dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dAlpha <= 0.0 || dBeta <= 0.0){
        return PENALTY;
    }
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        dLogP = logLikelihood(nA, dAlpha, dBeta);
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = logLikelihood(0, dAlpha, dBeta);
    
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
            
            if(verbose) printf("converged to minimum at\n");
        }
        
        if(verbose){
            printf ("%5d ", iter);
            
            for (i = 0; i < nP; i++) printf("%10.3e ", gsl_vector_get(s->x, i));
            
            printf("f() = %7.3f size = %.3f\n", s->fval, size);
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
void loadAbundance(t_Data *ptData, SAbundVector* rank)
{
    int nNA = 0; int  nL = 0, nJ = 0;
    int maxRank = rank->getMaxRank();
    
    for(int i = 1; i < maxRank; i++){ if (rank->get(i) != 0) { nNA++; } }

    int **aanAbund = NULL;
    aanAbund = (int **) malloc(nNA*sizeof(int*));
    
    int count = 0;
    for(int i = 1; i < maxRank; i++){
        
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
void freeAbundance(t_Data *ptData)
{
    for(int i = 0; i < ptData->nNA; i++){
        free(ptData->aanAbund[i]);
    }
    free(ptData->aanAbund);
}
/***********************************************************************/
void outputResults(gsl_vector *ptX, t_Data *ptData)
{
    double dAlpha = 0.0, dBeta = 0.0, dS = 0.0, dL = 0.0;
    
    dAlpha = gsl_vector_get(ptX, 0);
    
    dBeta  = gsl_vector_get(ptX, 1);
    
    dS = gsl_vector_get(ptX, 2);
    
    dL = nLogLikelihood(ptX, ptData);
    
    printf("\nML simplex: a = %.2f b = %.2f S = %.2f NLL = %.2f\n",dAlpha, dBeta, dS, dL);
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
    
    //printf("%e %e %e\n",dDeltaA,dDeltaB,dDeltaG);
    
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
    double dLogNot0 = 0.0, dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dAlpha <= 0.0 || dBeta <= 0.0){
        return PENALTY;
    }
    
    for(i = 0; i < ptData->nNA; i++){
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        dLogP = logLikelihood(nA, dAlpha, dBeta);
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = logLikelihood(0, dAlpha, dBeta);
    
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
    FILE    *sfp = NULL;
    int nS = 0, nSDash = 0, nIter = 0;
    double dRand = 0.0, dNLL = 0.0;
    void   *pvRet = NULL;
    
    /*set up random number generator*/
    T        = gsl_rng_default;
    ptGSLRNG = gsl_rng_alloc (T);
    
    nS = (int) floor(gsl_vector_get(ptX,2));
    
    dNLL = negLogLikelihood(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), nS,(void*) ptData);
    
    sprintf(szSampleFile,"%s_%d%s", ptParams->szOutFileStub.c_str(), ptMetroInit->nThread, SAMPLE_FILE_SUFFIX);
    
    sfp = fopen(szSampleFile, "w");
    if(!sfp){
        exit(EXIT_FAILURE);
    }
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        double dA = 0.0, dNLLDash = 0.0;
        
        getProposal(ptGSLRNG, ptXDash, ptX, &nSDash, nS, ptParams);
        
        dNLLDash = negLogLikelihood(gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), nSDash, (void*) ptData);
        //printf("X' %e %e %e %d %f\n", gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), gsl_vector_get(ptXDash,2), nSDash, dNLLDash);
        //printf("X %e %e %e %d %f\n", gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), gsl_vector_get(ptX,2), nS, dNLL);
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
            fprintf(sfp, "%d,%e,%e,%d,%f\n",nIter,gsl_vector_get(ptX, 0), gsl_vector_get(ptX, 1), nS, dNLL);
            fflush(sfp);
        }
        
        nIter++;
    }
    
    fclose(sfp);
    
    /*free up allocated memory*/
    gsl_vector_free(ptXDash);
    free(szSampleFile);
    gsl_rng_free(ptGSLRNG);
    
    return pvRet;
}

/***********************************************************************/
void writeThread(t_MetroInit *ptMetroInit)
{
    gsl_vector *ptX = ptMetroInit->ptX;
    printf("%d: a = %.2f b = %.2f S = %.2f\n", ptMetroInit->nThread,
           gsl_vector_get(ptX, 0),
           gsl_vector_get(ptX, 1),
           gsl_vector_get(ptX, 2));
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
    
    printf("\nMCMC iter = %d sigmaA = %.2f sigmaB = %.2f sigmaS = %.2f\n",
           ptParams->nIter, ptParams->dSigmaA, ptParams->dSigmaB, ptParams->dSigmaS);
    
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
    
    atMetroInit[1].ptParams = ptParams;
    atMetroInit[1].ptData   = ptData;
    atMetroInit[1].ptX      = ptX2;
    atMetroInit[1].nThread  = 1;
    atMetroInit[1].lSeed    = ptParams->lSeed + 1;
    atMetroInit[1].nAccepted = 0;
    
    atMetroInit[2].ptParams = ptParams;
    atMetroInit[2].ptData   = ptData;
    atMetroInit[2].ptX      = ptX3;
    atMetroInit[2].nThread  = 2;
    atMetroInit[2].lSeed    = ptParams->lSeed + 2;
    atMetroInit[2].nAccepted = 0;
    
    writeThread(&atMetroInit[0]);
    writeThread(&atMetroInit[1]);
    writeThread(&atMetroInit[2]);
    
    iret1 = pthread_create(&thread1, NULL, metropolis, (void*) &atMetroInit[0]);
    iret2 = pthread_create(&thread2, NULL, metropolis, (void*) &atMetroInit[1]);
    iret3 = pthread_create(&thread3, NULL, metropolis, (void*) &atMetroInit[2]);
    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);
    pthread_join(thread3, NULL);
    
    
    printf("%d: accept. ratio %d/%d = %f\n", atMetroInit[0].nThread,
           atMetroInit[0].nAccepted, ptParams->nIter,((double) atMetroInit[0].nAccepted)/((double) ptParams->nIter));
    
    printf("%d: accept. ratio %d/%d = %f\n", atMetroInit[1].nThread,
           atMetroInit[1].nAccepted, ptParams->nIter,((double) atMetroInit[1].nAccepted)/((double) ptParams->nIter));
    
    printf("%d: accept. ratio %d/%d = %f\n", atMetroInit[2].nThread,
           atMetroInit[2].nAccepted, ptParams->nIter, ((double) atMetroInit[2].nAccepted)/((double) ptParams->nIter));
    
    gsl_vector_free(ptX1); gsl_vector_free(ptX2); gsl_vector_free(ptX3);
}

#endif
/***********************************************************************/
double MetroIG::getValues(SAbundVector* rank){
    try {
        
        int  i = 0, nNA     = 0;
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaA = sigmaA; tParams.dSigmaB = sigmaB; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub;
        t_Data   tData;
#ifdef USE_GSL
        loadAbundance(&tData, rank);

        gsl_vector* ptX = gsl_vector_alloc(3); /*parameter estimates*/
        t_MetroInit atMetroInit[3];
        
        int maxRank = rank->getMaxRank();
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl

        gsl_rng_env_setup();
        
        gsl_set_error_handler_off();

        /*set initial estimates for parameters*/
        gsl_vector_set(ptX, 0, INIT_A);
        gsl_vector_set(ptX, 1, INIT_B);
        gsl_vector_set(ptX, 2, numOTUs*2);
        
        Chao1 chao; EstOutput results = chao.getValues(rank);
        
        printf("D = %d L = %d Chao = %f\n",numOTUs, sampled, results[0]);
        
        minimiseSimplex(ptX, 3, (void*) &tData, &nLogLikelihood);
        
        outputResults(ptX, &tData);
        
        if(tParams.nIter > 0){ mcmc(&tParams, &tData, ptX); }
        
        /*free up allocated memory*/
        gsl_vector_free(ptX);
      
        freeAbundance(&tData);
#endif    
        double result = 0;
        
        if (isnan(result) || isinf(result)) { result = 0; }
        
        return result;
    }
    catch(exception& e) {
        m->errorOut(e, "MetroIG", "getValues");
        exit(1);
    }
}
/***********************************************************************/
