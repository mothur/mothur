//
//  metrologstudent.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/2/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#include "metrologstudent.hpp"

/***********************************************************************/
MetroLogStudent::MetroLogStudent(int fi, double sigm, double sigv, double sign, double sigS, int n, string st) : sigmaM(sigm), sigmaV(sigv), sigmaN(sign), sigmaS(sigS), nIters(n), outFileStub(st), fitIters(fi), DiversityCalculator(false) {}
/***********************************************************************/


#ifdef USE_GSL
/***********************************************************************/
double nLogLikelihood2(const gsl_vector * x, void * params)
{
    double dMDash  = gsl_vector_get(x,0), dV = gsl_vector_get(x,1);
    double dNu  = gsl_vector_get(x,2);
    int    nS = (int) floor(gsl_vector_get(x, 3));
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dV <= 0.0 || dNu < 1.0){
        return PENALTY;
    }
    
    DiversityUtils dutils("metrols");
    
    for(i = 0; i < ptData->nNA; i++){
        
        if (dutils.m->getControl_pressed()) { break; }
        
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < 100){ //MAX_QUAD
            dLogP = dutils.logLikelihoodQuad(nA, dMDash, dV, dNu);
        }
        else{
            dLogP = dutils.logLikelihoodRampal(nA, dMDash, dV, dNu);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = dutils.logLikelihoodQuad(0, dMDash, dV, dNu);
    
    dLog1 = (nS - ptData->nL)*dLog0;
    
    dLog2 = - gsl_sf_lnfact(nS - ptData->nL);
    
    dLog3 = gsl_sf_lnfact(nS);
    
    dLogL += dLog1 + dLog2 + dLog3;
    
    /*return*/
    return -dLogL;
}
/***********************************************************************/
double negLogLikelihood(double dMDash, double dV, double dNu, int nS, void * params)
{
    t_Data *ptData = (t_Data *) params;
    int    i       = 0;
    double dLogL   = 0.0;
    double dLog0 = 0.0, dLog1 = 0.0, dLog2 = 0.0, dLog3 = 0.0;
    
    if(dV <= 0.0 || dNu < 1.0){
        return PENALTY;
    }
    DiversityUtils dutils("metrols");
    
    for(i = 0; i < ptData->nNA; i++){
        
        if (dutils.m->getControl_pressed()) { break; }
        
        double dLogP = 0.0;
        int    nA    = ptData->aanAbund[i][0];
        
        if(nA < 100){ //MAX_QUAD
            dLogP = dutils.logLikelihoodQuad(nA, dMDash, dV, dNu);
        }
        else{
            dLogP = dutils.logLikelihoodRampal(nA, dMDash, dV, dNu);
        }
        
        dLogL += ((double) ptData->aanAbund[i][1])*dLogP;
        
        dLogL -= gsl_sf_lnfact(ptData->aanAbund[i][1]);
        
    }
    
    dLog0 = dutils.logLikelihoodQuad(0, dMDash, dV, dNu);
    
    dLog1 = (nS - ptData->nL)*dLog0;
    
    dLog2 = - gsl_sf_lnfact(nS - ptData->nL);
    
    dLog3 = gsl_sf_lnfact(nS);
    
    dLogL += dLog1 + dLog2 + dLog3;
    
    /*return*/
    return -dLogL;
}
/***********************************************************************/
void* metropolis2 (void * pvInitMetro)
{
    t_MetroInit *ptMetroInit  = (t_MetroInit *) pvInitMetro;
    gsl_vector  *ptX          = ptMetroInit->ptX;
    t_Data      *ptData       = ptMetroInit->ptData;
    t_Params    *ptParams     = ptMetroInit->ptParams;
    gsl_vector  *ptXDash      = gsl_vector_alloc(4); /*proposal*/
    char *szSampleFile = (char *) malloc(1024*sizeof(char));
    const gsl_rng_type *T;
    gsl_rng            *ptGSLRNG;
    int nS = 0, nSDash = 0, nIter = 0;
    double dRand = 0.0, dNLL = 0.0;
    void   *pvRet = NULL;
    
    /*set up random number generator*/
    T        = gsl_rng_default;
    ptGSLRNG = gsl_rng_alloc (T);
    
    nS = (int) floor(gsl_vector_get(ptX,3));
    
    dNLL = negLogLikelihood(gsl_vector_get(ptX,0), gsl_vector_get(ptX,1), gsl_vector_get(ptX,2), nS,(void*) ptData);
    
    string filename = ptParams->szOutFileStub + "_" + toString(ptMetroInit->nThread) + ".sample";
    
    ofstream out; Utils util; util.openOutputFile(filename, out);
    out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
    
    /*seed random number generator*/
    gsl_rng_set(ptGSLRNG, ptMetroInit->lSeed);
    
    DiversityUtils dutils("metrols");
    
    /*now perform simple Metropolis algorithm*/
    while(nIter < ptParams->nIter){
        
        if (dutils.m->getControl_pressed()) { break; }
        
        double dA = 0.0, dNLLDash = 0.0;
        
        dutils.getProposal(ptGSLRNG, ptXDash, ptX, &nSDash, nS, ptParams);
        
        dNLLDash = negLogLikelihood(gsl_vector_get(ptXDash,0), gsl_vector_get(ptXDash,1), gsl_vector_get(ptXDash,2), nSDash, (void*) ptData);
        
        dA = exp(dNLL - dNLLDash);
        
        if(dA > 1.0){ dA = 1.0; }
        
        dRand = gsl_rng_uniform(ptGSLRNG);
        
        if(dRand < dA){
            gsl_vector_memcpy(ptX, ptXDash);
            nS = nSDash;
            dNLL = dNLLDash;
            ptMetroInit->nAccepted++;
        }
        
        if(nIter % 10 == 0){
            out << nIter << "," << gsl_vector_get(ptX, 0) << "," << gsl_vector_get(ptX, 1) << "," << gsl_vector_get(ptX, 2) << "," << nS << "," << dNLL << endl;
        }else if(nIter % 10000 == 0){
            cout << nIter << endl;
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
vector<string> MetroLogStudent::getValues(SAbundVector* rank){
    try {
        
        t_Params tParams; tParams.nIter = nIters; tParams.dSigmaX = sigmaM; tParams.dSigmaY = sigmaV; tParams.dSigmaN = sigmaN; tParams.dSigmaS = sigmaS; tParams.szOutFileStub = outFileStub; tParams.lSeed = m->getRandomSeed();
        t_Data   tData;
        
        int bestSample = 0;
        
#ifdef USE_GSL
        
        DiversityUtils dutils("metrols");
        
        dutils.loadAbundance(&tData, rank);
        
        int sampled = rank->getNumSeqs(); //nj
        int numOTUs = rank->getNumBins(); //nl
        
        gsl_vector* ptX = gsl_vector_alloc(4); /*parameter estimates*/
        
        gsl_rng_env_setup();
        
        gsl_set_error_handler_off();
        
        /*set initial estimates for parameters*/
        gsl_vector_set(ptX, 0, -10.0); //INIT_M
        gsl_vector_set(ptX, 1, 20.0); //INIT_V
        gsl_vector_set(ptX, 2, 20.0); //INIT_N
        gsl_vector_set(ptX, 3, numOTUs*2);
        
        double chaoResult = dutils.chao(&tData);
        m->mothurOut("\nMetroLogStudent - D = " + toString(numOTUs) + " L = " + toString(sampled) +  " Chao = " + toString(chaoResult) +  "\n");
        
        dutils.minimiseSimplex(ptX, 4, (void*) &tData, &nLogLikelihood2, 0.1, 1.0e-3, 100000);
        dutils.outputResults(ptX, &tData, &nLogLikelihood2);
        
        if(tParams.nIter > 0){
            
            vector<double> acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis2);
            
            if (fitIters != 0) {
                
                acceptRatioPos defaultRatio = dutils.findBest(acceptanceRates);
                
                int numTries = 1;
                map<double, acceptRatioPos> sigmaToAccept; //sigma value -> acceptance ratio
                map<acceptRatioPos, double> acceptToSigma; //acceptance ratio -> sigma value
                
                acceptRatioPos temp; //1.0 and pos 0 be default
                double newSigmaA = sigmaM/2.0;              //0.10
                sigmaToAccept[newSigmaA] = temp;            //0.05
                newSigmaA /= 2.0;
                sigmaToAccept[newSigmaA] = temp;            //0.025
                sigmaToAccept[(sigmaM/10.0)] = temp;        //0.01
                sigmaToAccept[(sigmaM/100.0)] = temp;       //0.001
                
                
                newSigmaA = sigmaM + (sigmaM/2.0);         //0.15
                sigmaToAccept[newSigmaA] = temp;
                newSigmaA = sigmaM + (sigmaM/10.0);         //0.11
                sigmaToAccept[newSigmaA] = temp;
                newSigmaA = sigmaM + (sigmaM/4.0);         //0.125
                sigmaToAccept[newSigmaA] = temp;
                newSigmaA = sigmaM + (3*sigmaM/4.0);        //0.175
                sigmaToAccept[newSigmaA] = temp;
                newSigmaA = sigmaM+sigmaM;                  //0.2
                sigmaToAccept[newSigmaA] = temp;

                
                for (map<double, acceptRatioPos>::iterator it = sigmaToAccept.begin(); it != sigmaToAccept.end(); it++) {
                    if (m->getControl_pressed()) { break; }
                    
                    tParams.dSigmaX = it->first; tParams.dSigmaY = it->first; tParams.dSigmaN = it->first;
                    
                    acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis2);
                    
                    it->second = dutils.findBest(acceptanceRates);
                    
                    acceptToSigma[it->second] = it->first;
                    
                    if (it->second.acceptRatio <= 0.05) { break; }
                }
                
                
                sigmaToAccept[sigmaM] = defaultRatio; //0.1
                acceptToSigma[defaultRatio] = sigmaM;
                
                //adjust around closest value
                acceptRatioPos thisBest = acceptToSigma.begin()->first;
                sigmaM = acceptToSigma.begin()->second;
                
                double factor = 0.05;
                
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
                    
                    acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis2);
                    
                    thisBest = dutils.findBest(acceptanceRates);
                    
                    acceptToSigma[thisBest] = tParams.dSigmaX;
                    sigmaToAccept[tParams.dSigmaX] = thisBest;
                    
                    numTries++;
                }
                
                if (numTries == fitIters) {
                    sigmaM = acceptToSigma.begin()->second;
                    tParams.dSigmaX = sigmaM; tParams.dSigmaY = sigmaM; tParams.dSigmaN = sigmaM;
                    
                    acceptanceRates = dutils.mcmc(&tParams, &tData, ptX, &metropolis2);
                    
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
        m->errorOut(e, "MetroLogStudent", "getValues");
        exit(1);
    }
}
/***********************************************************************/
