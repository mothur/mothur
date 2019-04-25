//
//  metroig.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef metroig_hpp
#define metroig_hpp

#include "mothurout.h"
#include "sabundvector.hpp"
#include "diversityutils.hpp"


typedef struct s_Params
{
    long lSeed;
    
    string szOutFileStub;
    
    double dSigmaA;
    
    double dSigmaB;
    
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
/***********************************************************************/

class MetroIG  {
    
public:
    MetroIG(double sigA, double sigB, double sigS, int n, string stub) : sigmaA(sigA), sigmaB(sigB), sigmaS(sigS), nIters(n), outFileStub(stub) { m = MothurOut::getInstance(); }
    
    vector<string> getValues(SAbundVector* rank);
    
    bool requiresSample() { return false; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    DiversityUtils dutils;
    
    double sigmaA, sigmaB, sigmaS;
    int nIters;
    string outFileStub;
    
    void loadAbundance(t_Data *ptData, SAbundVector* rank);
    void freeAbundance(t_Data *ptData);
    
    
};

/***********************************************************************/





#endif /* metroig_hpp */
