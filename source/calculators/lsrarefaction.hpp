//
//  lsrarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/20/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef lsrarefaction_hpp
#define lsrarefaction_hpp

#include "diversityutils.hpp"

/***********************************************************************/

class LSRarefaction   {
    
public:
    
    LSRarefaction(double c) : coverage(c) { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
    double coverage;
    
};

/***********************************************************************/




#endif /* lsrarefaction_hpp */
