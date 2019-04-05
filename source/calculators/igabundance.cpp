//
//  igabundance.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "igabundance.hpp"

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
#ifdef USE_GSL
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
#endif
/***********************************************************************/

double IGAbundance::getValues(SAbundVector* rank, int n) {
    try {
        double results = 0;
        
        if (isnan(results) || isinf(results)) { results= 0; }
        
        return results;
    }
    catch(exception& e) {
        m->errorOut(e, "IGAbundance", "getValues");
        exit(1);
    }
}
/***********************************************************************/


