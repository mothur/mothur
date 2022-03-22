//
//  metrosichel.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef metrosichel_hpp
#define metrosichel_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

//MetroSichel - Fits the compound Poisson Sichel dist
/***********************************************************************/

class MetroSichel : public DiversityCalculator {
    
public:
    
    MetroSichel(int af, double siga, double sigb, double sigg, double sigS, int n, string st);
    ~MetroSichel() = default;
    
    vector<string> getValues(SAbundVector* rank);
    
    string getTag() { return "si"; }
    
    
private:
    
    double sigmaA, sigmaB, sigmaG, sigmaS;
    int nIters, fitIters;
    string outFileStub;
    
};

/***********************************************************************/



#endif /* metrosichel_hpp */
