//
//  lsabundance.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/16/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef lsabundance_hpp
#define lsabundance_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class LSAbundance : public DiversityCalculator {
    
public:
    
    LSAbundance();
    
    vector<double> getValues(int mr, vector<mcmcSample>& sampling);
    
    string getTag() { return "ls"; }
    
private:
    
    
};

/***********************************************************************/



#endif /* lsabundance_hpp */
