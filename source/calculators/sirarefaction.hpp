//
//  sirarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef sirarefaction_hpp
#define sirarefaction_hpp

#include "diversityutils.hpp"

/***********************************************************************/

class SIRarefaction   {
    
public:
    
    SIRarefaction(double c) : coverage(c) { m = MothurOut::getInstance(); }
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
    
    double coverage;
    
};

/***********************************************************************/




#endif /* sirarefaction_hpp */
