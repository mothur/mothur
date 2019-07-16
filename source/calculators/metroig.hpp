//
//  metroig.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef metroig_hpp
#define metroig_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class MetroIG : public DiversityCalculator {
    
public:
    MetroIG(int fi, double sigA, double sigB, double sigS, int n, string stub);
    
    vector<string> getValues(SAbundVector* rank);
    
    string getTag() { return "ig"; }
    
    
private:
    
    double sigmaA, sigmaB, sigmaS;
    int nIters, fitIters;
    string outFileStub;
    
    
    
};

/***********************************************************************/

#endif /* metroig_hpp */


