//
//  kmerdist.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/29/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef kmerdist_hpp
#define kmerdist_hpp

#include "calculator.h"

/* This calculator is based on Edgar's method described here,
 Edgar, 2004
 Edgar, R. C. (2004).
 Muscle: a multiple sequence alignment method with reduced time and space complexity.
 BMC Bioinformatics, 5:113.
 */
/**************************************************************************************************/

class KmerDist {
    
public:
    
    KmerDist(int k); 
    
    double calcDist(Sequence A, Sequence B);
    string getCitation() { return "http://mothur.org"; }
    vector<double> calcDist(vector<kmerCount> A, vector<int> B, int); //A contains indexes to kmers it contains, B is size maxKmer intialized to false with kmers it contains set to true


private:
    int kmerSize, maxKmer;
    MothurOut* m;
};
/**************************************************************************************************/

#endif /* kmerdist_hpp */
