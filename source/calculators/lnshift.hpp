//
//  lnshift.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/14/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef lnshift_hpp
#define lnshift_hpp

#include "diversityutils.hpp"

/***********************************************************************/

class LNShift   {
    
public:
    
    LNShift() { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
};

/***********************************************************************/



#endif /* lnshift_hpp */
