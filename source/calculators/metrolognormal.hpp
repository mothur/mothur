//
//  metrolognormal.h
//  Mothur
//
//  Created by Sarah Westcott on 4/25/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef metrolognormal_h
#define metrolognormal_h

#include "diversityutils.hpp"


//MetroLogNormal - fits a compound Poisson Log-Normal distn to a sample
/***********************************************************************/

class MetroLogNormal  {
    
public:
    
    MetroLogNormal(double sigx, double sigy, double sigS, int n, string st) : sigmaX(sigx), sigmaY(sigy), sigmaS(sigS), nIters(n), outFileStub(st) { m = MothurOut::getInstance(); }
    ~MetroLogNormal() {}
    
    vector<string> getValues(SAbundVector* rank);
    
    bool requiresSample() { return false; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    
    double sigmaX, sigmaY, sigmaS;
    int nIters;
    string outFileStub;
    
};

/***********************************************************************/


#endif /* metrolognormal_h */
