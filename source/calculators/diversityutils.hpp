//
//  diversityutils.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversityutils_hpp
#define diversityutils_hpp

#include "mothurout.h"
#include "utils.hpp"

/***********************************************************************/

class DiversityUtils   {
    
public:
    DiversityUtils(){ m = MothurOut::getInstance(); }
    
    #ifdef USE_GSL
    
    double logLikelihood(int n, double dAlpha, double dBeta);
    int bessel(double* pdResult, int n, double dAlpha, double dBeta);
    double sd(int n, double dAlpha, double dBeta);
    
    #endif
    
    double f2X(double x, double dA, double dB, double dNDash);
    double fX(double x, double dA, double dB, double dNDash);
    
private:
    Utils util;
    MothurOut* m;
};

/***********************************************************************/



#endif /* diversityutils_hpp */
