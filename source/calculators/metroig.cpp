//
//  metroig.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metroig.hpp"

/*constants for simplex minimisation*/

/***********************************************************************/
MetroIG::MetroIG(int fi, double sigA, double sigB, double sigS, int n, string stub) : sigmaA(sigA), sigmaB(sigB), sigmaS(sigS), nIters(n), fitIters(fi), outFileStub(stub), DiversityCalculator(false) {}
/***********************************************************************/

#ifdef USE_GSL

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
    
    DiversityUtils dutils("metroig");
    
    for(i = 0; i < ptData->nNA; i++){
        
        if (dutils.m->getControl_pressed()) { break; }
        
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
    
    DiversityUtils dutils("metroig");
    
    for(i = 0; i < ptData->nNA; i++){
        
        if (dutils.m->getControl_pressed()) { break; }
        
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
    out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    DiversityUtils dutils("metroig");
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        double dA = 0.0, dNLLDash = 0.0;
        
        if (dutils.m->getControl_pressed()) { break; }
        
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
        
        DiversityUtils dutils("metroig");
        
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
        
        dutils.minimiseSimplex(ptX, 3, (void*) &tData, &nLogLikelihood0, 0.1, 1.0e-2, 100000);
        
        dutils.outputResults(ptX, &tData, &nLogLikelihood0);
        
        int bestSample = 0;
        
        if(tParams.nIter > 0){
    
            vector<double> acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis0);
            
            if (fitIters != 0) {
                
                acceptRatioPos defaultRatio = dutils.findBest(acceptanceRates);
                
                int numTries = 1;
                map<double, acceptRatioPos> sigmaToAccept; //sigma value -> acceptance ratio
                map<acceptRatioPos, double> acceptToSigma; //acceptance ratio -> sigma value
                
                acceptRatioPos temp; //1.0 and pos 0 be default
                sigmaToAccept[(sigmaA/10.0)] = temp; //0.01
                sigmaToAccept[(sigmaA/100.0)] = temp; //0.001
                
                double newSigmaA = sigmaA/2.0;
                sigmaToAccept[newSigmaA] = temp; //0.05
                newSigmaA /= 2.0;
                sigmaToAccept[newSigmaA] = temp; //0.025
                sigmaA /= 200.0;
                sigmaToAccept[newSigmaA] = temp; //0.0005
                
                for (map<double, acceptRatioPos>::iterator it = sigmaToAccept.begin(); it != sigmaToAccept.end(); it++) {
                    if (m->getControl_pressed()) { break; }
                    
                    tParams.dSigmaX = it->first; tParams.dSigmaY = it->first; tParams.dSigmaN = it->first;
                    
                    acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis0);
                    
                    it->second = dutils.findBest(acceptanceRates);
                    
                    acceptToSigma[it->second] = it->first;
                    
                    if (it->second.acceptRatio <= 0.05) { break; }
                }
                
                sigmaToAccept[sigmaA] = defaultRatio; //0.1
                acceptToSigma[defaultRatio] = sigmaA;
                
                //adjust around closest value
                acceptRatioPos thisBest = acceptToSigma.begin()->first;
                sigmaA = acceptToSigma.begin()->second;
                
                double factor = sigmaA / 2.0;
                
                while ((thisBest.acceptRatio > 0.05) && (numTries < fitIters)) {
                    if (m->getControl_pressed()) { break; }
                    
                    m->mothurOut("\nFit try: " + toString(numTries) + "\n");
                    
                    if (thisBest.acceptRatio < 0.45) {
                        
                        tParams.dSigmaX -= factor; tParams.dSigmaY -= factor; tParams.dSigmaN -= factor;
                        
                    }else if (thisBest.acceptRatio > 0.55) {
                        
                        tParams.dSigmaX += factor; tParams.dSigmaY += factor; tParams.dSigmaN += factor;
                        
                    }
                    
                    map<double, acceptRatioPos>::iterator it = sigmaToAccept.find(tParams.dSigmaX);
                    
                    if (it != sigmaToAccept.end()) { //we already tried this value, take average of 2 best tries
                        map<acceptRatioPos, double>::iterator it2 = acceptToSigma.begin();
                        double average = it2->second;
                        
                        it2++;
                        average += it2->second;
                        average /= 2.0;
                        
                        tParams.dSigmaX = average; tParams.dSigmaY = average; tParams.dSigmaN = average;
                    }
                    
                    acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis0);
                    
                    thisBest = dutils.findBest(acceptanceRates);
                    
                    acceptToSigma[thisBest] = tParams.dSigmaX;
                    sigmaToAccept[tParams.dSigmaX] = thisBest;
                    
                    numTries++;
                }
                
                if (numTries == fitIters) {
                    sigmaA = acceptToSigma.begin()->second;
                    tParams.dSigmaX = sigmaA; tParams.dSigmaY = sigmaA; tParams.dSigmaN = sigmaA;
                    
                    acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis0);
                    
                    thisBest = dutils.findBest(acceptanceRates);
                }
                
                bestSample = thisBest.pos;
            }

        }
        
        
        /*free up allocated memory*/
        gsl_vector_free(ptX);
        
        dutils.freeAbundance(&tData);
        
#endif
        
        outputs.push_back(outFileStub + "_" + toString(bestSample) + ".sample");
        if (bestSample == 0) {  outputs.push_back(outFileStub + "_1.sample"); outputs.push_back(outFileStub + "_2.sample");  }
        else if (bestSample == 1) {  outputs.push_back(outFileStub + "_0.sample"); outputs.push_back(outFileStub + "_2.sample");  }
        else if (bestSample == 2) {  outputs.push_back(outFileStub + "_0.sample"); outputs.push_back(outFileStub + "_1.sample");  }
        
        return outputs;
    }
    catch(exception& e) {
        m->errorOut(e, "MetroIG", "getValues");
        exit(1);
    }
}
/***********************************************************************/

