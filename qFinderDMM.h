//
//  qFinderDMM.h
//  pds_dmm
//
//  Created by Patrick Schloss on 11/8/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pds_dmm_qFinderDMM_h
#define pds_dmm_qFinderDMM_h

#include "communitytype.h"

/**************************************************************************************************/

class qFinderDMM : public CommunityTypeFinder {
  
public:
    qFinderDMM(vector<vector<int> >, int);
    void printFitData(ofstream&);
    void printFitData(ostream&, double);
    
private:
   
    void optimizeLambda();
    void calculatePiK();

    double negativeLogEvidenceLambdaPi(vector<double>&);
    void negativeLogDerivEvidenceLambdaPi(vector<double>&, vector<double>&);
    double getNegativeLogEvidence(vector<double>&, int);
    double getNegativeLogLikelihood();
    
    
    int lineMinimizeFletcher(vector<double>&, vector<double>&, double, double, double, double&, double&, vector<double>&, vector<double>&);
    int bfgs2_Solver(vector<double>&);//, double, double);
    
   

        
};

/**************************************************************************************************/

#endif
