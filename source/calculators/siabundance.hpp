//
//  siabundance.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/22/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef siabundance_hpp
#define siabundance_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class SIAbundance : public DiversityCalculator  {
    
public:
    
    SIAbundance();
    
    vector<double> getValues(int mr, vector<mcmcSample>& sampling);
    
    string getTag() { return "si"; }
    
private:
    
    
};

/***********************************************************************/



#endif /* siabundance_hpp */
