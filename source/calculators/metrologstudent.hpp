//
//  metrologstudent.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/2/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef metrologstudent_hpp
#define metrologstudent_hpp

#include "diversityutils.hpp"


//MetroLogStudent - Samples the Poisson Log-Student distn to SADs
/***********************************************************************/

class MetroLogStudent  {
    
public:
    
    MetroLogStudent(double sigm, double sigv, double sign, double sigS, int n, string st) : sigmaM(sigm), sigmaV(sigv), sigmaN(sign), sigmaS(sigS), nIters(n), outFileStub(st) { m = MothurOut::getInstance(); }
    ~MetroLogStudent() {}
    
    vector<string> getValues(SAbundVector* rank);
    
    bool requiresSample() { return false; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    
    double sigmaM, sigmaV, sigmaN, sigmaS;
    int nIters;
    string outFileStub;
    
};

/***********************************************************************/


#endif /* metrologstudent_hpp */
