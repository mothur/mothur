//
//  diversityutils.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "diversityutils.hpp"

/***********************************************************************/
void DiversityUtils::loadAbundance(t_Data *ptData, SAbundVector* rank)
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
double derivExponent(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
    double dTemp = (x - dMDash)/dV, dRet = 0.0;
    
    dRet = ((double) n) - exp(x) - dTemp;
    
    return dRet;
}

/***********************************************************************/
double DiversityUtils::chao(t_Data *ptData)
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
void DiversityUtils::freeAbundance(t_Data *ptData)
{
    for(int i = 0; i < ptData->nNA; i++){
        free(ptData->aanAbund[i]);
    }
    free(ptData->aanAbund);
}
/***********************************************************************/
double DiversityUtils::fX(double x, double dA, double dB, double dNDash)
{
    double dTemp1 = (dA*(x - dB)*(x - dB))/x;
    
    return log(x) - (1.0/dNDash)*(x + dTemp1);
}
/***********************************************************************/
double DiversityUtils::f2X(double x, double dA, double dB, double dNDash)
{
    double dRet = 0.0, dTemp = 2.0*dA*dB*dB;
    
    dRet = (1.0/(x*x))*(1.0 + (1.0/dNDash)*(dTemp/x));
    
    return -dRet;
}
 #ifdef USE_GSL

/***********************************************************************/
double fMu(double x, void* pvParams)
{
    t_IGParams* ptIGParams = (t_IGParams *) pvParams;
    
    DiversityUtils dutils("igrarefaction");

    double dAlphaDD = ptIGParams->dAlpha*sqrt(x);
    double dBetaDD  = ptIGParams->dBeta*x;
    double dLogP0   = dutils.logLikelihood(0, dAlphaDD, dBetaDD);
    
    return (1.0 - exp(dLogP0)) - ptIGParams->dC;
}

/***********************************************************************/
int DiversityUtils::solveF(double x_lo, double x_hi, void* params, double tol, double *xsolve)
{
    int status, iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    gsl_function F;
    
    F.function = derivExponent;
    if (method == "igrarefaction") { F.function = fMu;  }
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
double DiversityUtils::sd(int n, double dAlpha, double dBeta)
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
int DiversityUtils::bessel(double* pdResult, int n, double dAlpha, double dBeta)
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
double DiversityUtils::logLikelihood(int n, double dAlpha, double dBeta)
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
double DiversityUtils::sd(int n, double dAlpha, double dBeta, double dGamma)
{
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
int DiversityUtils::bessel(double* pdResult, int n, double dAlpha, double dBeta, double dGamma)
{
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
double DiversityUtils::logLikelihood(int n, double dAlpha, double dBeta, double dGamma)
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
    
    status = bessel(&dRet,n, dAlpha,dBeta,dGamma);
    if(status == FALSE){
        dRet = sd(n, dAlpha,dBeta,dGamma);
    }
    
    return dRet - dLogFacN;
}

/***********************************************************************/

void DiversityUtils::outputResults(gsl_vector *ptX, t_Data *ptData, double (*f)(const gsl_vector*, void* params))
{
    double dAlpha = 0.0, dBeta = 0.0, dS = 0.0, dL = 0.0, dGamma = 0.0;
    
    dAlpha = gsl_vector_get(ptX, 0);
    
    dBeta  = gsl_vector_get(ptX, 1);
    
    if ((method == "metrols") || (method == "metrosichel"))  {
        dGamma = gsl_vector_get(ptX, 2);
        dS = gsl_vector_get(ptX, 3);
    }else {
        dS = gsl_vector_get(ptX, 2);
    }
    
    dL = f(ptX, ptData);
    
    if (method == "metroig")        { m->mothurOut("\nMetroIG - ML simplex: a = " + toString(dAlpha) +  " b = " + toString(dBeta) +  " S = " + toString(dS) +  " NLL = " + toString(dL) + "\n");            }
    else if (method == "metroln")   { m->mothurOut("\nMetroLogNormal - ML simplex: M = " + toString(dAlpha) +  " V = " + toString(dBeta) +  " S = " + toString(dS) +  " NLL = " + toString(dL) + "\n");     }
    else if (method == "metrols")   { m->mothurOut("\nMetroLogStudent - ML simplex: M = " + toString(dAlpha) +  " V = " + toString(dBeta) +  " Nu = " + toString(dGamma) +  " S = " + toString(dS) +  " NLL = " + toString(dL) + "\n");     }
    else if (method == "metrosichel")   { m->mothurOut("\nMetroSichel - ML simplex: a = " + toString(dAlpha) +  " b = " + toString(dBeta) +  " g = " + toString(dGamma) +  " S = " + toString(dS) +  " NLL = " + toString(dL) + "\n");     }
    
}
/***********************************************************************/
int DiversityUtils::minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params), double initSimplexSize, double minSimplexSize, double maxSimplexSize)
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
    
    if (method == "metroig")        {
        gsl_vector_set_all(ss, initSimplexSize);
        gsl_vector_set(ss,nP - 1,0.1*gsl_vector_get(ptX,0));
    }
    else if (method == "metroln")   {
        gsl_vector_set_all(ss, initSimplexSize);
        gsl_vector_set(ss,2,0.1*gsl_vector_get(ptX,2));
    }
    else if (method == "metrols" ) {
        for(i = 0; i < nP; i++){
            gsl_vector_set(ss, i,initSimplexSize*fabs(gsl_vector_get(ptX,i)));
        }
    }else if (method == "metrosichel" ) {
        for(i = 0; i < nP; i++){
            gsl_vector_set(ss, i,initSimplexSize*gsl_vector_get(ptX,i));
        }
    }
    
    /* Initialize method and iterate */
    minex_func.f = f;
    minex_func.n = nP;
    minex_func.params = pvData;
    
    s = gsl_multimin_fminimizer_alloc (T, nP);
    gsl_multimin_fminimizer_set(s, &minex_func, ptX, ss);
    
    do{
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if(status) { break; }
        
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, minSimplexSize);
        
        if(status == GSL_SUCCESS){
            for(i = 0; i < nP; i++){
                gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
            }
        }
    }
    while(status == GSL_CONTINUE && iter < maxSimplexSize);
    
    if(status == GSL_CONTINUE){
        for(i = 0; i < nP; i++){
            gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
        }
    }
    
    
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    return status;
}
/***********************************************************************/
void DiversityUtils::getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams)
{
    double dDeltaS =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaS);
    double dDeltaA =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaX);
    double dDeltaB =  gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaY);
    
    double dDeltaN =  0;
    int    nSDash = 0;
    if ((method == "metrols") || (method == "metrosichel")) { dDeltaN = gsl_ran_gaussian(ptGSLRNG, ptParams->dSigmaN); }
    
    gsl_vector_set(ptXDash, 0, gsl_vector_get(ptX,0) + dDeltaA);
    gsl_vector_set(ptXDash, 1, gsl_vector_get(ptX,1) + dDeltaB);
    if ((method == "metrols") || (method == "metrosichel")) {  gsl_vector_set(ptXDash, 2, gsl_vector_get(ptX,2) + dDeltaN); }
    
    nSDash = nS + (int) floor(dDeltaS);
    if(nSDash < 1){
        nSDash = 1;
    }
    (*pnSDash) = nSDash;
}

/***********************************************************************/
void DiversityUtils::mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX, void* f (void * pvInitMetro))
{
    int ptXSize = 3;
    if ((method == "metrols") || (method == "metrosichel")) {  ptXSize = 4;  }
    
    pthread_t thread1, thread2, thread3;
    int       iret1  , iret2  , iret3;
    gsl_vector *ptX1 = gsl_vector_alloc(ptXSize),
    *ptX2 = gsl_vector_alloc(ptXSize),
    *ptX3 = gsl_vector_alloc(ptXSize);
    t_MetroInit atMetroInit[3];
    
    if (method == "metrols") {  m->mothurOut("\nMCMC iter = " + toString(ptParams->nIter) + " sigmaM = " + toString(ptParams->dSigmaX) +  " sigmaV = " + toString(ptParams->dSigmaY) +  " sigmaN = " + toString(ptParams->dSigmaN) +  " sigmaS = " + toString(ptParams->dSigmaS) + "\n"); }
    else if (method == "metrosichel") {  m->mothurOut("\nMCMC iter = " + toString(ptParams->nIter) + " sigmaA = " + toString(ptParams->dSigmaX) +  " sigmaB = " + toString(ptParams->dSigmaY) +  " sigmaG = " + toString(ptParams->dSigmaN) +  " sigmaS = " + toString(ptParams->dSigmaS) + "\n"); }
    else { m->mothurOut("\nMCMC iter = " + toString(ptParams->nIter) + " sigmaX = " + toString(ptParams->dSigmaX) +  " sigmaY = " + toString(ptParams->dSigmaY) +  " sigmaS = " + toString(ptParams->dSigmaS) + "\n"); }
    
    gsl_vector_memcpy(ptX1, ptX);
    
    gsl_vector_set(ptX2, 0, gsl_vector_get(ptX,0) + 2.0*ptParams->dSigmaX);
    gsl_vector_set(ptX2, 1, gsl_vector_get(ptX,1) + 2.0*ptParams->dSigmaY);
    
    if ((method == "metrols") || (method == "metrosichel")) {
        gsl_vector_set(ptX2, 2, gsl_vector_get(ptX,2) + 2.0*ptParams->dSigmaN);
        gsl_vector_set(ptX2, 3, gsl_vector_get(ptX,3) + 2.0*ptParams->dSigmaS);
    }
    else { gsl_vector_set(ptX2, 2, gsl_vector_get(ptX,2) + 2.0*ptParams->dSigmaS);  }
    
    
    gsl_vector_set(ptX3, 0, gsl_vector_get(ptX,0) - 2.0*ptParams->dSigmaX);
    gsl_vector_set(ptX3, 1, gsl_vector_get(ptX,1) - 2.0*ptParams->dSigmaY);
    
    if ((method == "metrols") || (method == "metrosichel")) {
        gsl_vector_set(ptX3, 2, gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaN);
        
        if(gsl_vector_get(ptX,3) - 2.0*ptParams->dSigmaS > (double) ptData->nL){
            gsl_vector_set(ptX3, 3, gsl_vector_get(ptX,3) - 2.0*ptParams->dSigmaS);   }
        else{ gsl_vector_set(ptX3, 3, (double) ptData->nL);   }
    }else {
        
        if(gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaS > (double) ptData->nL){ gsl_vector_set(ptX3, 2, gsl_vector_get(ptX,2) - 2.0*ptParams->dSigmaS); }
        else{ gsl_vector_set(ptX3, 2, (double) ptData->nL); }
    }
    
    atMetroInit[0].ptParams = ptParams;
    atMetroInit[0].ptData   = ptData;
    atMetroInit[0].ptX      = ptX1;
    atMetroInit[0].nThread  = 0;
    atMetroInit[0].lSeed    = ptParams->lSeed;
    atMetroInit[0].nAccepted = 0;
    
    //write thread 0
    
    if ((method == "metrols") || (method == "metrosichel")) { m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX1, 0)) +  " b = " + toString(gsl_vector_get(ptX1, 1)) +  " g = " + toString(gsl_vector_get(ptX1, 2)) +  " S = " + toString(gsl_vector_get(ptX1, 3)) + "\n"); }
    else { m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX1, 0)) +  " b = " + toString(gsl_vector_get(ptX1, 1)) +  " S = " + toString(gsl_vector_get(ptX1, 2)) + "\n"); }
    
    atMetroInit[1].ptParams = ptParams;
    atMetroInit[1].ptData   = ptData;
    atMetroInit[1].ptX      = ptX2;
    atMetroInit[1].nThread  = 1;
    atMetroInit[1].lSeed    = ptParams->lSeed + 1;
    atMetroInit[1].nAccepted = 0;
    
    //write thread 1
    if ((method == "metrols") || (method == "metrosichel")) { m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX2, 0)) +  " b = " + toString(gsl_vector_get(ptX2, 1)) +  " g = " + toString(gsl_vector_get(ptX2, 2)) +  " S = " + toString(gsl_vector_get(ptX2, 3)) + "\n"); }
    else { m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX2, 0)) +  " b = " + toString(gsl_vector_get(ptX2, 1)) +  " S = " + toString(gsl_vector_get(ptX2, 2)) + "\n"); }
    
    atMetroInit[2].ptParams = ptParams;
    atMetroInit[2].ptData   = ptData;
    atMetroInit[2].ptX      = ptX3;
    atMetroInit[2].nThread  = 2;
    atMetroInit[2].lSeed    = ptParams->lSeed + 2;
    atMetroInit[2].nAccepted = 0;
    
    //write thread 2
    if ((method == "metrols") || (method == "metrosichel")) { m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX3, 0)) +  " b = " + toString(gsl_vector_get(ptX3, 1)) +  " g = " + toString(gsl_vector_get(ptX3, 2)) +  " S = " + toString(gsl_vector_get(ptX3, 3)) + "\n"); }
    else { m->mothurOut(toString(atMetroInit[0].nThread) + ": a = " + toString(gsl_vector_get(ptX3, 0)) +  " b = " + toString(gsl_vector_get(ptX3, 1)) +  " S = " + toString(gsl_vector_get(ptX3, 2)) + "\n"); }
    
    iret1 = pthread_create(&thread1, NULL, f, (void*) &atMetroInit[0]);
    iret2 = pthread_create(&thread2, NULL, f, (void*) &atMetroInit[1]);
    iret3 = pthread_create(&thread3, NULL, f, (void*) &atMetroInit[2]);
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






