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
#include "diversitycalc.h"

//MetroLogNormal - fits a compound Poisson Log-Normal distn to a sample
/***********************************************************************/

class MetroLogNormal : public DiversityCalculator {
    
public:
    
    MetroLogNormal(int fi, double sigx, double sigy, double sigS, int n, string st);
    ~MetroLogNormal() {}
    
    vector<string> getValues(SAbundVector* rank);
    
    string getTag() { return "ln"; }
    
private:
    
    double sigmaX, sigmaY, sigmaS;
    int nIters, fitIters;
    string outFileStub;
    
};

/***********************************************************************/


#endif /* metrolognormal_h */
