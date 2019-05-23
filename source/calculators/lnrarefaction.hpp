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
#include "diversitycalc.h"

/***********************************************************************/

class LNRarefaction : public DiversityCalculator   {
    
public:
    
    LNRarefaction(double c);// : coverage(c) { m = MothurOut::getInstance(); }
    
    vector<double> getValues(int ns, vector<mcmcSample>& sampling);
    
    
private:
    
    double coverage;
    
};

/***********************************************************************/



#endif /* lnrarefaction_hpp */
