//
//  igrarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/6/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef igrarefaction_hpp
#define igrarefaction_hpp

#include "diversityutils.hpp"


//IGRarefaction
/***********************************************************************/

class IGRarefaction  {
    
public:
    
    IGRarefaction(double c) : coverage(c) { m = MothurOut::getInstance(); }
    ~IGRarefaction() {}
    
    vector<double> getValues(SAbundVector* rank, vector<mcmcSample>& sampling);
    
    bool requiresSample() { return true; }
    
    
private:
    
    Utils util;
    MothurOut* m;
    
    double coverage;
    
    double calcMu(t_IGParams *ptIGParams);
    
    
};

/***********************************************************************/



#endif /* igrarefaction_hpp */
