//
//  metroig.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metroig.hpp"

/*constants for simplex minimisation*/
//#define PENALTY           1.0e20

//#define INIT_A            1.0
//#define INIT_B            5.0
//#define INIT_SIMPLEX_SIZE 0.1
//#define INIT_S_SS         0.1
//#define MIN_SIMPLEX_SIZE  1.0e-2
//#define MAX_SIMPLEX_ITER  100000
//#define MAX_LINE_LENGTH   1024
//#define SAMPLE_FILE_SUFFIX ".sample"

//#define DEF_SIGMA     0.1
//#define DEF_SIGMA_S   100.0

//#define DEF_ITER   250000
//#define SLICE      10
//#define PRECISION  1.0e-10

#ifdef USE_GSL


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

#endif
/***********************************************************************/
vector<string> MetroIG::getValues(SAbundVector* rank){
    try {
    
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaX = sigmaA; tParams.dSigmaY = sigmaB; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub; tParams.lSeed = m->getRandomSeed();
        t_Data   tData;
#ifdef USE_GSL
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
        m->mothurOut("\nD = " + toString(numOTUs) + " L = " + toString(sampled) +  " Chao = " + toString(chaoResult) +  "\n");
        
        minimiseSimplex(ptX, 3, (void*) &tData, &nLogLikelihood);
        
        outputResults(ptX, &tData);
        
        if(tParams.nIter > 0){
           dutils.mcmc(&tParams, &tData, ptX);
        }
        
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
