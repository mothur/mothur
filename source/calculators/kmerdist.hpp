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

/**************************************************************************************************/

class KmerDist {
    
public:
    
    KmerDist(int k); 
    
    double calcDist(Sequence A, Sequence B);
    double calcDist(vector<int> A, vector<bool> B, int); //A contains indexes to kmers it contains, B is size maxKmer intialized to false with kmers it contains set to true

private:
    int kmerSize, maxKmer;
    MothurOut* m;
};
/**************************************************************************************************/

#endif /* kmerdist_hpp */
