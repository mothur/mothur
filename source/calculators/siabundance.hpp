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

/***********************************************************************/

class SIAbundance   {
    
public:
    
    SIAbundance() { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
};

/***********************************************************************/



#endif /* siabundance_hpp */
