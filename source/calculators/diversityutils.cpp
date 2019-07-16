//
//  diversityutils.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "diversityutils.hpp"
/***********************************************************************/
double f1_2(double x, void *pvParams)
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
double f1(double x, void *pvParams)
{
    t_LNParams *ptLNParams = (t_LNParams *) pvParams;
    double dMDash = ptLNParams->dMDash, dV = ptLNParams->dV, n = ptLNParams->n;
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
void DiversityUtils::loadAbundance(t_Data *ptData, SAbundVector* rank) {
    try {
        int nNA = 0; int  nL = 0, nJ = 0;
        int maxRank = rank->getMaxRank();
        
        for(int i = 1; i <= maxRank; i++){ if (rank->get(i) != 0) { nNA++; } }
        
        int **aanAbund = NULL;
        aanAbund = (int **) malloc(nNA*sizeof(int*));
        
        int count = 0;
        for(int i = 1; i <= maxRank; i++){
            
            if (m->getControl_pressed()) { break; }
            
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
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "loadAbundance");
        exit(1);
    }
}
/***********************************************************************/
void DiversityUtils::freeAbundance(t_Data *ptData) {
    try {
        for(int i = 0; i < ptData->nNA; i++){ free(ptData->aanAbund[i]); }
        free(ptData->aanAbund);
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "freeAbundance");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::chao(t_Data *ptData) {
    try {
        double n1 = 0.0, n2 = 0.0;
        int **aanAbund = ptData->aanAbund;
        
        if(aanAbund[0][0] == 1 && aanAbund[1][0] == 2){
            
            n1 = (double) aanAbund[0][1];
            n2 = (double) aanAbund[1][1];
            
            return ((double) ptData->nL) + 0.5*((n1*n1)/n2);
        }
        else{ return -1.0; }
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "chao");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::fX(double x, double dA, double dB, double dNDash){
    try {
        double dTemp1 = (dA*(x - dB)*(x - dB))/x;
        
        return log(x) - (1.0/dNDash)*(x + dTemp1);
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "fX");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::f2X(double x, double dA, double dB, double dNDash) {
    try {
        double dRet = 0.0, dTemp = 2.0*dA*dB*dB;
        
        dRet = (1.0/(x*x))*(1.0 + (1.0/dNDash)*(dTemp/x));
        
        return -dRet;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "f2X");
        exit(1);
    }
}
 #ifdef USE_GSL

/***********************************************************************/
double fMu_sirarefaction(double x, void* pvParams)
{
    DiversityUtils dutils("sirarefaction");
    
    t_LSParams* ptSIParams = (t_LSParams *) pvParams;
    double dAlphaDD = ptSIParams->dMDash*sqrt(x);
    double dBetaDD  = ptSIParams->dV*x;
    double dLogP0   = dutils.logLikelihood(0, dAlphaDD, dBetaDD, ptSIParams->dNu);
    
    return (1.0 - exp(dLogP0)) - ptSIParams->dC;
}
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
double fMu_lsrarefaction(double x, void* pvParams)
{
    
    DiversityUtils dutils("lsrarefaction");
    
    t_LSParams* ptLSParams = (t_LSParams *) pvParams;
    double dMDD            = ptLSParams->dMDash + x;
    double dLogP0          = dutils.logLikelihoodQuad(0, dMDD, ptLSParams->dV, ptLSParams->dNu);
    
    return (1.0 - exp(dLogP0)) - ptLSParams->dC;
}
/***********************************************************************/
double fMu_lnrarefaction(double x, void* pvParams)
{
    t_IGParams* ptIGParams = (t_IGParams *) pvParams;
    
    DiversityUtils dutils("lnrarefaction");
    
    double dMDD = ptIGParams->dAlpha + x;
    double dLogP0 = dutils.logLikelihoodQuad(0, dMDD, ptIGParams->dBeta);
    
    return (1.0 - exp(dLogP0)) - ptIGParams->dC;
}
/***********************************************************************/
double DiversityUtils::calcMu(void *pvParams){
    try {
        double dLogMu = 0.0;
        
        if (method == "lnrarefaction") {
            t_IGParams* ptIGParams = (t_IGParams *) pvParams;
            solveF(0, 1.0e7, ptIGParams, 1.0e-7, &dLogMu);
            return exp(dLogMu);
        }else if ((method == "igrarefaction") || (method == "sirarefaction")) {
            t_IGParams* ptIGParams = (t_IGParams *) pvParams;
            solveF(1.0, 1.0e10, ptIGParams, 1.0e-7, &dLogMu);
            return dLogMu;
        }else if (method == "lsrarefaction") {
            t_LSParams *ptLSParams = (t_LSParams *) pvParams;
            solveF(0, 1.0e7, ptLSParams, 1.0e-7, &dLogMu);
            return exp(dLogMu);
        }
        
        return dLogMu;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "calcMu");
        exit(1);
    }
}
//***********************************************************************/
double DiversityUtils::logStirlingsGamma(double dZ) { return 0.5*log(2.0*M_PI) + (dZ - 0.5)*log(dZ) - dZ; }
/***********************************************************************/
double DiversityUtils::logLikelihoodQuad(int n, double dMDash, double dV){
    try {
        gsl_integration_workspace *ptGSLWS =
        gsl_integration_workspace_alloc(1000);
        double dLogFac1   = 0.0, dLogFacN  = 0.0;
        double dResult = 0.0, dError = 0.0, dPrecision = 0.0;
        gsl_function tGSLF;
        double dEst = dMDash + ((double)n)*dV, dA = 0.0, dB = 0.0;
        
        t_LNParams tLNParams; tLNParams.n = n; tLNParams.dMDash = dMDash; tLNParams.dV = dV;
        
        tGSLF.function = &f1;
        if (method == "metrols") {  tGSLF.function = &f1Log;  }
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
            
            if ((method == "metroln") || (method == "metrols")) { if(fabs(dUpper) > 1.0e-7){  solveF(0.0, dUpper, (void *) &tLNParams, 1.0e-5, &dMax); } }
            else {  solveF(0.0, dUpper, derivExponent, (void *) &tLNParams, 1.0e-5, &dMax);  } //lnabund, lnshift, lnrarefact
            
            dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
            
            dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
        }
        else{
            double dMax = 0.0;
            double dLower = dEst - dV;
            double dUpper = (((double) n) + (dMDash/dV) - 1.0)/(1.0 + 1.0/dV);
            double dVar   = 0.0;
            
            if ((method == "metroln") || (method == "metrols")) {
                if(fabs(dUpper - dLower) > 1.0e-7){ solveF(dLower, dUpper, (void *) &tLNParams, 1.0e-5, &dMax); }
                else{ dMax = 0.5*(dLower + dUpper); }
            }else {  solveF(dLower, dUpper, derivExponent, (void *) &tLNParams, 1.0e-5, &dMax); } //lnabund, lnshift, lnrarefact
            
            dVar = sqrt(1.0/((1.0/dV) + exp(dMax)));
            
            dA = dMax - V_MULT*dVar; dB = dMax + V_MULT*dVar;
        }
        
        if(n < 10)  {  dPrecision = HI_PRECISION;   }
        else        {  dPrecision = LO_PRECISION;   }
        
        gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
        
        gsl_integration_workspace_free(ptGSLWS);
        
        return log(dResult) - dLogFacN -0.5*dLogFac1;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "logLikelihoodQuad");
        exit(1);
    }
}
//***********************************************************************/
double DiversityUtils::logLikelihoodQuad(int n, double dMDash, double dV, double dNu){
    try {
        gsl_integration_workspace *ptGSLWS =
        gsl_integration_workspace_alloc(1000);
        double dLogFac1   = 0.0, dLogFacN  = 0.0;
        double dN = (double) n, dResult = 0.0, dError = 0.0, dPrecision = 0.0;
        gsl_function tGSLF;
        t_LSParams tLSParams;
        double dA = 0.0, dB = 0.0;
        
        tLSParams.n = n; tLSParams.dMDash = dMDash; tLSParams.dV = dV; tLSParams.dNu = dNu;
        
        tGSLF.function = &f1_2;
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
        
        if(n < 10)  { dPrecision = HI_PRECISION; }
        else        { dPrecision = LO_PRECISION; }
        
        gsl_integration_qag(&tGSLF, dA, dB, dPrecision, 0.0, 1000, GSL_INTEG_GAUSS61, ptGSLWS, &dResult, &dError);
        
        //printf("%f %f\n", dResult, dError);
        
        gsl_integration_workspace_free(ptGSLWS);
        
        return log(dResult) - dLogFacN + dLogFac1 - 0.5*log(dV);
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "logLikelihoodQuad");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::logLikelihoodRampal(int n, double dMDash, double dV){
    try {
        double dN = (double) n;
        double dLogLik = 0.0, dTemp = gsl_pow_int(log(dN) - dMDash,2), dTemp3 = gsl_pow_int(log(dN) - dMDash,3);
        
        dLogLik = -0.5*log(2.0*M_PI*dV) - log(dN) - (dTemp/(2.0*dV));
        
        dLogLik += log(1.0 + 1.0/(2.0*dN*dV)*(dTemp/dV + log(dN) - dMDash - 1.0)
                       + 1.0/(6.0*dN*dN*dV*dV*dV)*(3.0*dV*dV - (3.0*dV - 2.0*dV*dV)*(dMDash - log(dN))
                                                   - 3.0*dV*dTemp + dTemp3));
        
        return dLogLik;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "logLikelihoodRampal");
        exit(1);
    }
}
//***********************************************************************/
double DiversityUtils::logLikelihoodRampal(int n, double dMDash, double dV, double dNu){
    try {
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
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "logLikelihoodRampal");
        exit(1);
    }
}
//***********************************************************************/
int DiversityUtils::solveF(double x_lo, double x_hi, double (*f)(double, void*), void* params, double tol, double *xsolve){
    try {
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
            
            if (m->getControl_pressed()) { break; }
            
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
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "solveF");
        exit(1);
    }
}
//***********************************************************************/
int DiversityUtils::solveF(double x_lo, double x_hi, void* params, double tol, double *xsolve){
    try {
        int status, iter = 0, max_iter = 100;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        double r = 0;
        gsl_function F;
        
        F.function = &derivExponent;
        if (method == "igrarefaction") {  F.function = &fMu_igrarefaction; }
        else if (method == "lnrarefaction") {  F.function = &fMu_lnrarefaction; }
        else if (method == "lsrarefaction") {  F.function = &fMu_lsrarefaction; }
        else if (method == "sirarefaction") {  F.function = &fMu_sirarefaction; }
        F.params = params;
        
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        gsl_root_fsolver_set (s, &F, x_lo, x_hi);
        
        do{
            iter++;
            
            if (m->getControl_pressed()) { break; }
            
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
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "solveF");
        exit(1);
    }
}

/***********************************************************************/
double DiversityUtils::logLikelihood(int n, double dAlpha, double dBeta){
    try {
        double dLogFacN = 0.0;
        bool status      = false;
        double dRet     = 0.0;
        
        if(n < 50){
            dLogFacN = gsl_sf_fact(n);
            dLogFacN = log(dLogFacN);
        }
        else{
            dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
        }
        
        status = bessel(&dRet,n, dAlpha,dBeta);
        
        if(!status){ dRet = sd(n, dAlpha,dBeta); }
        
        return dRet - dLogFacN;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "logLikelihood");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::logLikelihood(int n, double dAlpha, double dBeta, double dGamma){
    try {
        double dLogFacN = 0.0;
        bool status      = false;
        double dRet     = 0.0;
        
        if(n < 50){
            dLogFacN = gsl_sf_fact(n);
            dLogFacN = log(dLogFacN);
        }
        else{
            dLogFacN = gsl_sf_lngamma(((double) n) + 1.0);
        }
        
        status = bessel(&dRet,n, dAlpha,dBeta,dGamma);
        
        if(!status){ dRet = sd(n, dAlpha,dBeta,dGamma); }
        
        return dRet - dLogFacN;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "logLikelihood");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::sd(int n, double dAlpha, double dBeta){
    try {
        double dGamma = -0.5;
        double dA = 0.5*(-1.0 + sqrt(1.0 + (dAlpha*dAlpha)/(dBeta*dBeta)));
        double dN = (double) n, dNDash = dN + dGamma - 1.0, dRN = 1.0/dN;
        double dTemp1 = (0.5*dN)/(1.0 + dA), dTemp2 = 4.0*dRN*dRN*(1.0 + dA)*dA*dBeta*dBeta;
        double dXStar = dTemp1*(1.0 + sqrt(1.0 + dTemp2));
        double dFX = fX(dXStar, dA, dBeta, dNDash);
        double d2FX = -dNDash*f2X(dXStar, dA, dBeta, dNDash);
        double dLogK = 0.0, dGamma1 = dGamma;
        
        if(dGamma1 < 0.0){ dGamma1 *= -1.0; } //invert sign
        
        dLogK = gsl_sf_bessel_lnKnu(dGamma1,2.0*dA*dBeta);
        
        return -2.0*dA*dBeta -log(2.0) -dLogK -dGamma*log(dBeta) + dNDash*dFX + 0.5*log(2.0*M_PI) - 0.5*log(d2FX);
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "sd");
        exit(1);
    }
}
/***********************************************************************/
double DiversityUtils::sd(int n, double dAlpha, double dBeta, double dGamma){
    try {
        double dA = 0.5*(-1.0 + sqrt(1.0 + (dAlpha*dAlpha)/(dBeta*dBeta)));
        double dN = (double) n, dNDash = dN + dGamma - 1.0, dRN = 1.0/dN;
        double dTemp1 = (0.5*dN)/(1.0 + dA), dTemp2 = 4.0*dRN*dRN*(1.0 + dA)*dA*dBeta*dBeta;
        double dXStar = dTemp1*(1.0 + sqrt(1.0 + dTemp2));
        double dFX = fX(dXStar, dA, dBeta, dNDash);
        double d2FX = -dNDash*f2X(dXStar, dA, dBeta, dNDash);
        double dLogK = 0.0, dGamma1 = dGamma;
        
        if(dGamma1 < 0.0){ dGamma1 *= -1.0; } //invert sign
        
        dLogK = gsl_sf_bessel_lnKnu(dGamma1,2.0*dA*dBeta);
        
        return -2.0*dA*dBeta -log(2.0) -dLogK -dGamma*log(dBeta) + dNDash*dFX + 0.5*log(2.0*M_PI) - 0.5*log(d2FX);
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "sd");
        exit(1);
    }
}
/***********************************************************************/
bool DiversityUtils::bessel(double* pdResult, int n, double dAlpha, double dBeta){
    try {
        double dGamma  = -0.5;
        double dResult = 0.0;
        double dOmega = 0.0, dGamma2 = 0.0;
        double dLogK1 = 0.0, dLogK2 = 0.0;
        double dN = (double) n, dNu = dGamma + dN;
        double dTemp1 = 0.0;
        
        if(dNu < 0.0){ dNu = -dNu; }
        
        if(dGamma < 0.0)    { dGamma2 = -dGamma;    }
        else                { dGamma2 = dGamma;     }
        
        dOmega = sqrt(dBeta*dBeta + dAlpha*dAlpha) - dBeta;
        
        dLogK2 = gsl_sf_bessel_lnKnu(dNu, dAlpha);
        
        if(!gsl_finite(dLogK2)){
            if(dAlpha < 0.1*sqrt(dNu + 1.0)){
                dLogK2 = gsl_sf_lngamma(dNu) + (dNu - 1.0)*log(2.0) - dNu*log(dAlpha);
            }else{
                (*pdResult) = dResult;
                return false;
            }
        }
        
        dLogK1 = dGamma*log(dOmega/dAlpha) -gsl_sf_bessel_lnKnu(dGamma2,dOmega);
        
        dTemp1 = log((dBeta*dOmega)/dAlpha);
        
        dResult = dN*dTemp1 + dLogK2 + dLogK1;
        
        (*pdResult) = dResult;
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "bessel");
        exit(1);
    }
}
/***********************************************************************/
bool DiversityUtils::bessel(double* pdResult, int n, double dAlpha, double dBeta, double dGamma){
    try {
        double dResult = 0.0;
        double dOmega = 0.0, dGamma2 = 0.0;
        double dLogK1 = 0.0, dLogK2 = 0.0;
        double dN = (double) n, dNu = dGamma + dN;
        double dTemp1 = 0.0;
        
        if(dNu < 0.0){ dNu = -dNu; }
        
        if(dGamma < 0.0){ dGamma2 = -dGamma;    }
        else            { dGamma2 = dGamma;     }
        
        dOmega = sqrt(dBeta*dBeta + dAlpha*dAlpha) - dBeta;
        
        dLogK2 = gsl_sf_bessel_lnKnu(dNu, dAlpha);
        
        if(!gsl_finite(dLogK2)){
            if(dAlpha < 0.1*sqrt(dNu + 1.0)){
                dLogK2 = gsl_sf_lngamma(dNu) + (dNu - 1.0)*log(2.0) - dNu*log(dAlpha);
            }else{
                (*pdResult) = dResult;
                return false;
            }
        }
        
        dLogK1 = dGamma*log(dOmega/dAlpha) -gsl_sf_bessel_lnKnu(dGamma2,dOmega);
        
        dTemp1 = log((dBeta*dOmega)/dAlpha);
        
        dResult = dN*dTemp1 + dLogK2 + dLogK1;
        
        (*pdResult) = dResult;
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "bessel");
        exit(1);
    }
}

/***********************************************************************/

void DiversityUtils::outputResults(gsl_vector *ptX, t_Data *ptData, double (*f)(const gsl_vector*, void* params)){
    try {
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
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "outputResults");
        exit(1);
    }
}
/***********************************************************************/
int DiversityUtils::minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params), double initSimplexSize, double minSimplexSize, double maxSimplexSize){
    try {
        const gsl_multimin_fminimizer_type *T =
        gsl_multimin_fminimizer_nmsimplex;
        gsl_multimin_fminimizer *s = NULL;
        gsl_vector *ss;
        gsl_multimin_function minex_func;
        size_t iter = 0;
        int status;
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
            for(int i = 0; i < nP; i++){
                gsl_vector_set(ss, i,initSimplexSize*fabs(gsl_vector_get(ptX,i)));
            }
        }else if (method == "metrosichel" ) {
            for(int i = 0; i < nP; i++){
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
            
            if (m->getControl_pressed()) { break; }
            
            status = gsl_multimin_fminimizer_iterate(s);
            
            if(status) { break; }
            
            size = gsl_multimin_fminimizer_size(s);
            status = gsl_multimin_test_size(size, minSimplexSize);
            
            if(status == GSL_SUCCESS){
                for(int i = 0; i < nP; i++){
                    gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
                }
            }
        }
        while(status == GSL_CONTINUE && iter < maxSimplexSize);
        
        if(status == GSL_CONTINUE){
            for(int i = 0; i < nP; i++){
                gsl_vector_set(ptX, i, gsl_vector_get(s->x, i));
            }
        }
        
        
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free (s);
        
        return status;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "minimiseSimplex");
        exit(1);
    }
}
/***********************************************************************/
void DiversityUtils::getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams){
    try {
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
        
        if(nSDash < 1){ nSDash = 1; }
        
        (*pnSDash) = nSDash;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "getProposal");
        exit(1);
    }
}
/***********************************************************************/
int DiversityUtils::fitSigma(vector<double> acceptanceRates, double sigmaA, int fitIters, t_Params *ptParams, t_Data *ptData, gsl_vector* ptX, void* f (void * pvInitMetro)){
    try {
        acceptRatioPos defaultRatio = findBest(acceptanceRates);
        
        if (defaultRatio.acceptRatio <= 0.05) { return defaultRatio.pos;  }
        
        int numTries = 1;
        map<double, acceptRatioPos> sigmaToAccept; //sigma value -> acceptance ratio
        map<acceptRatioPos, double> acceptToSigma; //acceptance ratio -> sigma value
        
        acceptRatioPos temp; //1.0 and pos 0 be default
        sigmaToAccept[(sigmaA/10.0)] = temp;        //0.01
        sigmaToAccept[(sigmaA/100.0)] = temp;       //0.001
        sigmaToAccept[(sigmaA/1000.0)] = temp;       //0.0001
        sigmaToAccept[(sigmaA/10000.0)] = temp;       //0.00001
        
        double newSigmaA = sigmaA + (sigmaA/2.0);         //0.15
        sigmaToAccept[newSigmaA] = temp;
        newSigmaA = sigmaA+sigmaA;                  //0.2
        sigmaToAccept[newSigmaA] = temp;
    
        //adjust around closest "high" and closest "low" values
        acceptRatioPos thisBestHigh, thisBestLow;
        map<acceptRatioPos, double> acceptToSigmaHigh; //acceptance ratio -> sigma value
        map<acceptRatioPos, double> acceptToSigmaLow; //acceptance ratio -> sigma value
        
        //set iters to 1000, get close to value then run with nIters
        int savedIters = ptParams->nIter;
        ptParams->nIter = 1000;
    
        for (map<double, acceptRatioPos>::iterator it = sigmaToAccept.begin(); it != sigmaToAccept.end(); it++) {
            if (m->getControl_pressed()) { break; }
            
            ptParams->dSigmaX = it->first; ptParams->dSigmaY = it->first; ptParams->dSigmaN = it->first;
            
            acceptanceRates = mcmc(ptParams, ptData, ptX, f);
            
            it->second = findBest(acceptanceRates);
            
            if (it->second.high) { //high
                if (it->second.acceptRatio < thisBestHigh.acceptRatio) {  thisBestHigh = it->second; }
                acceptToSigmaHigh[it->second] = it->first;
            }else { //low
                if (it->second.acceptRatio < thisBestLow.acceptRatio) {  thisBestLow = it->second; }
                acceptToSigmaLow[it->second] = it->first;
            }
            acceptToSigma[it->second] = it->first;
            
            if (it->second.acceptRatio <= 0.05) {
                //try with nIters to confirm
                ptParams->nIter = savedIters;
                
                acceptanceRates = mcmc(ptParams, ptData, ptX, f);
                
                it->second = findBest(acceptanceRates);
                
                if (it->second.high) { //high
                    if (it->second.acceptRatio < thisBestHigh.acceptRatio) {  thisBestHigh = it->second; }
                    acceptToSigmaHigh[it->second] = it->first;
                }else { //low
                    if (it->second.acceptRatio < thisBestLow.acceptRatio) {  thisBestLow = it->second; }
                    acceptToSigmaLow[it->second] = it->first;
                }
                acceptToSigma[it->second] = it->first;
                
                //if good value
                if (it->second.acceptRatio <= 0.05) { return it->second.pos;  }
                else {  ptParams->nIter = 1000;  }
            }
        }
        
        sigmaToAccept[sigmaA] = defaultRatio; //0.1
        acceptToSigma[defaultRatio] = sigmaA;
        
        double factor = 0.0; bool badHigh = false; bool badLow = false; double badFactor = 0.0;
        
        //find best high and check
        map<acceptRatioPos, double>::iterator itFind = acceptToSigma.find(thisBestHigh);
        if (itFind != acceptToSigma.end()) {
            if (thisBestHigh.acceptRatio > 0.25) {
                badHigh = true; badFactor += itFind->second;
            }else {
                factor += itFind->second;
                sigmaA = itFind->second;
            }
        }//else no high values
        
        //find best low and check
        itFind = acceptToSigma.find(thisBestLow);
        if (itFind != acceptToSigma.end()) {
            if (thisBestLow.acceptRatio > 0.25) { //below 25% acceptance, lets disregard
                badLow = true; badFactor += itFind->second;
            }else {
                factor += itFind->second;
                if (badHigh) { sigmaA = itFind->second; }
                else { if (sigmaA > itFind->second) { sigmaA = itFind->second; } }
            }
        }//no low values
        
        if (badHigh && badLow) {
            double increment = badFactor / (double)(fitIters);
            sigmaA = acceptToSigma.begin()->second; //sigma for best try
            sigmaA -= (increment*(fitIters/(double)2.0));
            factor = badFactor / (double)(fitIters);
        }else if (badHigh || badLow)  {
            double increment = factor / (double)(fitIters);
            sigmaA -= (increment*(fitIters/(double)2.0));
        }else { //good high and low
            double increment = factor / (double)(fitIters);
            sigmaA -= (increment*(fitIters/(double)2.0));
            factor /= (double) fitIters;
        }

        ptParams->dSigmaX = sigmaA; ptParams->dSigmaY = sigmaA; ptParams->dSigmaN = sigmaA;
        ptParams->nIter = savedIters;
        
        while ((thisBestLow.acceptRatio > 0.05) && (numTries < fitIters)) {
            if (m->getControl_pressed()) { break; }
            
            m->mothurOut("\nFit try: " + toString(numTries) + "\n");
            
            ptParams->dSigmaX += factor; ptParams->dSigmaY += factor; ptParams->dSigmaN += factor;
            
            map<double, acceptRatioPos>::iterator it = sigmaToAccept.find(ptParams->dSigmaX);
            
            if (it == sigmaToAccept.end()) {
                acceptanceRates = mcmc(ptParams, ptData, ptX, f);
            
                thisBestLow = findBest(acceptanceRates);
            
                acceptToSigma[thisBestLow] = ptParams->dSigmaX;
                sigmaToAccept[ptParams->dSigmaX] = thisBestLow;
                numTries++;
            }
        }
        
        if (numTries == fitIters) {
            sigmaA = acceptToSigma.begin()->second;
            ptParams->dSigmaX = sigmaA; ptParams->dSigmaY = sigmaA; ptParams->dSigmaN = sigmaA;
            
            acceptanceRates = mcmc(ptParams, ptData, ptX, f);
            
            thisBestLow = findBest(acceptanceRates);
        }
        
        if ((thisBestLow.acceptRatio > 0.05)) { m->mothurOut("\n[ERROR]: Unable to reach acceptable ratio, please review and set sigma parameters manually.\n"); m->setControl_pressed(true); }
        
        return thisBestLow.pos;
        
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "fitSigma");
        exit(1);
    }
}
/***********************************************************************/
vector<double> DiversityUtils::mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX, void* f (void * pvInitMetro)){
    try {
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
        if ((method == "metrols") || (method == "metrosichel")) { m->mothurOut(toString(atMetroInit[1].nThread) + ": a = " + toString(gsl_vector_get(ptX2, 0)) +  " b = " + toString(gsl_vector_get(ptX2, 1)) +  " g = " + toString(gsl_vector_get(ptX2, 2)) +  " S = " + toString(gsl_vector_get(ptX2, 3)) + "\n"); }
        else { m->mothurOut(toString(atMetroInit[1].nThread) + ": a = " + toString(gsl_vector_get(ptX2, 0)) +  " b = " + toString(gsl_vector_get(ptX2, 1)) +  " S = " + toString(gsl_vector_get(ptX2, 2)) + "\n"); }
        
        atMetroInit[2].ptParams = ptParams;
        atMetroInit[2].ptData   = ptData;
        atMetroInit[2].ptX      = ptX3;
        atMetroInit[2].nThread  = 2;
        atMetroInit[2].lSeed    = ptParams->lSeed + 2;
        atMetroInit[2].nAccepted = 0;
        
        //write thread 2
        if ((method == "metrols") || (method == "metrosichel")) { m->mothurOut(toString(atMetroInit[2].nThread) + ": a = " + toString(gsl_vector_get(ptX3, 0)) +  " b = " + toString(gsl_vector_get(ptX3, 1)) +  " g = " + toString(gsl_vector_get(ptX3, 2)) +  " S = " + toString(gsl_vector_get(ptX3, 3)) + "\n"); }
        else { m->mothurOut(toString(atMetroInit[2].nThread) + ": a = " + toString(gsl_vector_get(ptX3, 0)) +  " b = " + toString(gsl_vector_get(ptX3, 1)) +  " S = " + toString(gsl_vector_get(ptX3, 2)) + "\n"); }
        
        iret1 = pthread_create(&thread1, NULL, f, (void*) &atMetroInit[0]);
        iret2 = pthread_create(&thread2, NULL, f, (void*) &atMetroInit[1]);
        iret3 = pthread_create(&thread3, NULL, f, (void*) &atMetroInit[2]);
        pthread_join(thread1, NULL);
        pthread_join(thread2, NULL);
        pthread_join(thread3, NULL);
        
        m->mothurOut(toString(atMetroInit[0].nThread) +": accept. ratio " + toString(atMetroInit[0].nAccepted) + "/" + toString(ptParams->nIter) +  " = " + toString(((double) atMetroInit[0].nAccepted)/((double) ptParams->nIter)) +  "\n");
        m->mothurOut(toString(atMetroInit[1].nThread) +": accept. ratio " + toString(atMetroInit[1].nAccepted) + "/" + toString(ptParams->nIter) +  " = " + toString(((double) atMetroInit[1].nAccepted)/((double) ptParams->nIter)) +  "\n");
        m->mothurOut(toString(atMetroInit[2].nThread) +": accept. ratio " + toString(atMetroInit[2].nAccepted) + "/" + toString(ptParams->nIter) +  " = " + toString(((double) atMetroInit[2].nAccepted)/((double) ptParams->nIter)) +  "\n");
        
        vector<double> results;
        results.push_back(atMetroInit[0].nAccepted/((double) ptParams->nIter));
        results.push_back(atMetroInit[1].nAccepted/((double) ptParams->nIter));
        results.push_back(atMetroInit[2].nAccepted/((double) ptParams->nIter));
        
        gsl_vector_free(ptX1); gsl_vector_free(ptX2); gsl_vector_free(ptX3);
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "mcmc");
        exit(1);
    }
}

#endif

/***********************************************************************/
acceptRatioPos DiversityUtils::findBest(vector<double> acceptanceRates){
    try {
        
        double defaultSigmaAcc = fabs(0.5 - acceptanceRates[0]);  //"0" version
        bool high = true;
        if ((0.5 - acceptanceRates[0]) > 0.0) { high = false; }
        acceptRatioPos defaultRatio(defaultSigmaAcc, 0, high);
        
        
        if (defaultRatio.acceptRatio > fabs(0.5 - acceptanceRates[1])) {  //is the "1" version better?
            defaultRatio.acceptRatio = fabs(0.5 - acceptanceRates[1]);
            defaultRatio.pos = 1;
            defaultRatio.high = true;
            if ((0.5 - acceptanceRates[1]) > 0.0) { defaultRatio.high = false; }
        }
        if (defaultRatio.acceptRatio > fabs(0.5 - acceptanceRates[2])) {  //is the "2" version better?
            defaultRatio.acceptRatio = fabs(0.5 - acceptanceRates[2]);
            defaultRatio.pos = 2;
            defaultRatio.high = true;
            if ((0.5 - acceptanceRates[2]) > 0.0) { defaultRatio.high = false; }
        }
        
        return defaultRatio;
    }
    catch(exception& e) {
        m->errorOut(e, "DiversityUtils", "findBest");
        exit(1);
    }
}

/***********************************************************************/






