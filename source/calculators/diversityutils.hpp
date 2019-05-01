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

#define MIN_SIMPLEX_SIZE  1.0e-2
#define MAX_SIMPLEX_ITER  100000

/***********************************************************************/

class DiversityUtils   {
    
public:
    DiversityUtils(){ m = MothurOut::getInstance(); }
    
    #ifdef USE_GSL
    
    double logLikelihood(int n, double dAlpha, double dBeta);
    int bessel(double* pdResult, int n, double dAlpha, double dBeta);
    double sd(int n, double dAlpha, double dBeta);
    int minimiseSimplex(gsl_vector* ptX, size_t nP, void* pvData, double (*f)(const gsl_vector*, void* params), double, string);
    void mcmc(t_Params *ptParams, t_Data *ptData, gsl_vector* ptX, void* f (void * pvInitMetro));
    
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
