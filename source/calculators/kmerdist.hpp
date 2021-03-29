//
//  kmerdist.hpp
//  Mothur
//
//  Created by John Westcott on 3/29/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#ifndef kmerdist_hpp
#define kmerdist_hpp

#include "calculator.h"

/**************************************************************************************************/

class KmerDist : public DistCalc {
    
public:
    
    KmerDist(double c, int k); 
    
    double calcDist(Sequence A, Sequence B);

private:
    int kmerSize, maxKmer;
};
/**************************************************************************************************/

#endif /* kmerdist_hpp */
