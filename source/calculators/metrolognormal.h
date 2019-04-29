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
#include "diversityutils.hpp"

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
