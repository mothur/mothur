//
//  diversityutils.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversityutils_hpp
#define diversityutils_hpp

#define HI_PRECISION 1.0e-12
#define LO_PRECISION 1.0e-7

#define V_MULT          25.0
#define PENALTY           1.0e20
#define SLICE      10


#include "mothurout.h"
#include "utils.hpp"
#include "sabundvector.hpp"

/***********************************************************************/

class DiversityUtils   {
    
public:
    DiversityUtils(string met){ m = MothurOut::getInstance(); method = met; }
    
    #ifdef USE_GSL
    
    double logLikelihood(int n, double dAlpha, double dBeta);
    int bessel(double* pdResult, int n, double dAlpha, double dBeta);
    double sd(int n, double dAlpha, double dBeta);
    double logLikelihood(int n, double dAlpha, double dBeta, double);
    int bessel(double* pdResult, int n, double dAlpha, double dBeta, double);
    double sd(int n, double dAlpha, double dBeta, double);
    int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params), double, double, double);
    void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX, void* f (void * pvInitMetro));
    void outputResults(gsl_vector *ptX, t_Data *ptData, double (*f)(const gsl_vector*, void* params));
    void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams);
    int solveF(double x_lo, double x_hi, void* params, double tol, double *xsolve);
    double logLikelihoodRampal(int n, double dMDash, double dV);
    double logLikelihoodQuad(int n, double dMDash, double dV);
    
    
    #endif
    
    double f2X(double x, double dA, double dB, double dNDash);
    double fX(double x, double dA, double dB, double dNDash);
    double chao(t_Data *ptData);
    
    void loadAbundance(t_Data *ptData, SAbundVector* rank);
    void freeAbundance(t_Data *ptData);
    
    
private:
    
    Utils util;
    MothurOut* m;
    string method;
    
};

/***********************************************************************/



#endif /* diversityutils_hpp */
