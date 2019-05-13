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


//MetroSichel - Fits the compound Poisson Sichel dist
/***********************************************************************/

class MetroSichel  {
    
public:
    
    MetroSichel(double siga, double sigb, double sigg, double sigS, int n, string st) : sigmaA(siga), sigmaB(sigb), sigmaG(sigg), sigmaS(sigS), nIters(n), outFileStub(st) { m = MothurOut::getInstance(); }
    ~MetroSichel() {}
    
    vector<string> getValues(SAbundVector* rank);
    
    bool requiresSample() { return false; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    
    double sigmaA, sigmaB, sigmaG, sigmaS;
    int nIters;
    string outFileStub;
    
};

/***********************************************************************/



#endif /* metrosichel_hpp */
