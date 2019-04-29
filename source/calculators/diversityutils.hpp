//
//  diversityutils.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversityutils_hpp
#define diversityutils_hpp

#include "mothurout.h"
#include "utils.hpp"
#include "sabundvector.hpp"


typedef struct s_Params
{
    long lSeed;
    
    //char *szInputFile;
    
    string szOutFileStub;
    
    double dSigmaX;
    
    double dSigmaY;
    
    double dSigmaS;
    
    int nIter;
} t_Params;

typedef struct s_Data
{
    int nNA;
    
    int **aanAbund;
    
    int nL;
    
    int nJ;
}t_Data;

#ifdef USE_GSL
typedef struct s_MetroInit
{
    t_Params *ptParams;
    
    t_Data   *ptData;
    
    gsl_vector* ptX;
    
    int nAccepted;
    
    long lSeed;
    
    int nThread;
    
} t_MetroInit;
#endif

#ifdef USE_GSL

    double nLogLikelihood(const gsl_vector * x, void * params);
    void* metropolis (void * pvInitMetro);
    int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params)); //,  nLogLikelihood

#endif
/***********************************************************************/

class DiversityUtils   {
    
public:
    DiversityUtils(){ m = MothurOut::getInstance(); }
    
    #ifdef USE_GSL
    
    double logLikelihood(int n, double dAlpha, double dBeta);
    int bessel(double* pdResult, int n, double dAlpha, double dBeta);
    double sd(int n, double dAlpha, double dBeta);
    
    double negLogLikelihood(double dAlpha, double dBeta, int nS, void * params);
    void getProposal(gsl_rng *ptGSLRNG, gsl_vector *ptXDash, gsl_vector *ptX, int* pnSDash, int nS, t_Params *ptParams);
    void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX);
    
    #endif
    
    double f2X(double x, double dA, double dB, double dNDash);
    double fX(double x, double dA, double dB, double dNDash);
    double chao(t_Data *ptData);
    
    void loadAbundance(t_Data *ptData, SAbundVector* rank);
    void freeAbundance(t_Data *ptData);
    
private:
    Utils util;
    MothurOut* m;
};

/***********************************************************************/



#endif /* diversityutils_hpp */
