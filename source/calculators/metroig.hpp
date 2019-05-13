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

/***********************************************************************/

class MetroIG  {
    
public:
    MetroIG(double sigA, double sigB, double sigS, int n, string stub) : sigmaA(sigA), sigmaB(sigB), sigmaS(sigS), nIters(n), outFileStub(stub) { m = MothurOut::getInstance(); }
    
    vector<string> getValues(SAbundVector* rank);
    
    bool requiresSample() { return false; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    
    double sigmaA, sigmaB, sigmaS;
    int nIters;
    string outFileStub;
    
};

/***********************************************************************/

#endif /* metroig_hpp */


