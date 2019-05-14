//
//  lnrarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/13/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef lnrarefaction_hpp
#define lnrarefaction_hpp

#include "diversityutils.hpp"

/***********************************************************************/

class LNRarefaction   {
    
public:
    
    LNRarefaction(double c) : coverage(c) { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
    double coverage;
    
};

/***********************************************************************/



#endif /* lnrarefaction_hpp */
