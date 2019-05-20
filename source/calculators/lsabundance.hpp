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

/***********************************************************************/

class LSAbundance   {
    
public:
    
    LSAbundance() { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
};

/***********************************************************************/



#endif /* lsabundance_hpp */
