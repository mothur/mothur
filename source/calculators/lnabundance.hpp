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

/***********************************************************************/

class LNAbundance   {
    
public:
    
    LNAbundance() { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
};

/***********************************************************************/



#endif /* lnabundace_hpp */
