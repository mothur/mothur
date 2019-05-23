//
//  lnabundace.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/8/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef lnabundace_hpp
#define lnabundace_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class LNAbundance : public DiversityCalculator  {
    
public:
    
    LNAbundance();
    
    vector<double> getValues(int mr, vector<mcmcSample>& sampling);
    
    
private:
    
    
    
};

/***********************************************************************/



#endif /* lnabundace_hpp */
