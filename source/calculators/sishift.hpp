//
//  sishift.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef sishift_hpp
#define sishift_hpp

#include "diversityutils.hpp"

/***********************************************************************/

class SIShift   {
    
public:
    
    SIShift() { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
};

/***********************************************************************/




#endif /* sishift_hpp */
