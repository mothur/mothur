//
//  metrolognormal.h
//  Mothur
//
//  Created by Sarah Westcott on 4/25/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef metrolognormal_h
#define metrolognormal_h

#include "mothurout.h"

typedef struct s_Params
{
    long lSeed;
    
    char *szInputFile;
    
    char *szOutFileStub;
    
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

typedef struct s_MetroInit
{
    t_Params *ptParams;
    
    t_Data   *ptData;
    
    gsl_vector* ptX;
    
    int nAccepted;
    
    long lSeed;
    
    int nThread;
    
} t_MetroInit;


//MetroLogNormal - fits a compound Poisson Log-Normal distn to a sample
/***********************************************************************/

class MetroLogNormal  {
    
public:
    MetroLogNormal(double sigx, double sigy, double sigS, int n, string stub) : sigmaX(sigx), sigmaY(sigy), sigmaS(sigS), nIters(n), outFileStub(stub) { m = MothurOut::getInstance(); }
    
    vector<string> getValues(SAbundVector* rank);
    
    bool requiresSample() { return false; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    DiversityUtils dutils;
    
    double sigmaX, sigmaY, sigmaS;
    int nIters;
    string outFileStub;
    
    void loadAbundance(t_Data *ptData, SAbundVector* rank);
    void freeAbundance(t_Data *ptData);
    
    
};

/***********************************************************************/


#endif /* metrolognormal_h */
