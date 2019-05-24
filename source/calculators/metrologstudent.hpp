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
#include "diversitycalc.h"

//MetroLogStudent - Samples the Poisson Log-Student distn to SADs
/***********************************************************************/

class MetroLogStudent : public DiversityCalculator {
    
public:
    
    MetroLogStudent(double sigm, double sigv, double sign, double sigS, int n, string st);
    ~MetroLogStudent() {}
    
    vector<string> getValues(SAbundVector* rank);
    
    string getTag() { return "ls"; }
    
private:
    
    double sigmaM, sigmaV, sigmaN, sigmaS;
    int nIters;
    string outFileStub;
    
};

/***********************************************************************/


#endif /* metrologstudent_hpp */
