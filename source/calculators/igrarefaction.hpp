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
#include "diversitycalc.h"

//IGRarefaction
/***********************************************************************/

class IGRarefaction : public DiversityCalculator {
    
public:
    
    IGRarefaction(double c); 
    ~IGRarefaction() = default;
    
    vector<double> getValues(int ns, vector<mcmcSample>& sampling);
    
    string getTag() { return "ig"; }
    
private:
    
    double coverage;
    
    
};

/***********************************************************************/



#endif /* igrarefaction_hpp */
