//
//  qFinderDMM.h
//  pds_dmm
//
//  Created by Patrick Schloss on 11/8/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pds_dmm_qFinderDMM_h
#define pds_dmm_qFinderDMM_h

/**************************************************************************************************/

#include "mothurout.h"

/**************************************************************************************************/

class qFinderDMM {
  
public:
    qFinderDMM(vector<vector<int> >, int);
    double getNLL()     {    return currNLL;        }  
    double getAIC()     {    return aic;            }
    double getBIC()     {    return bic;            }
    double getLogDet()  {    return logDeterminant; }
    double getLaplace() {    return laplace;        }
    void printZMatrix(string, vector<string>);
    void printRelAbund(string, vector<string>);

private:
    MothurOut* m;
    void kMeans();
    void optimizeLambda();
    void calculatePiK();

    double negativeLogEvidenceLambdaPi(vector<double>&);
    void negativeLogDerivEvidenceLambdaPi(vector<double>&, vector<double>&);
    double getNegativeLogEvidence(vector<double>&, int);
    double getNegativeLogLikelihood();
    vector<vector<double> > getHessian();
    
    int lineMinimizeFletcher(vector<double>&, vector<double>&, double, double, double, double&, double&, vector<double>&, vector<double>&);
    int bfgs2_Solver(vector<double>&);//, double, double);
    double cheb_eval(const double[], int, double);
    double psi(double);
    double psi1(double);

    vector<vector<int> > countMatrix;
    vector<vector<double> > zMatrix;
    vector<vector<double> > lambdaMatrix;
    vector<double> weights;
    vector<vector<double> > error;
    
    int numPartitions;
    int numSamples;
    int numOTUs;
    int currentPartition;
    
    double currNLL;
    double aic;
    double bic;
    double logDeterminant;
    double laplace;
    
};

/**************************************************************************************************/

#endif
