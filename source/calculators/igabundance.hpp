//
//  igabundance.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef igabundance_hpp
#define igabundance_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class IGAbundance : public DiversityCalculator {
    
public:
    
    IGAbundance();
    
    vector<double> getValues(int mr, vector<mcmcSample>& sampling);
    
    string getTag() { return "ig"; }
    
private:
    Utils util;
    MothurOut* m;
    
};

/***********************************************************************/



#endif /* igabundance_hpp */
