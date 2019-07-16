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
#include "diversitycalc.h"

/***********************************************************************/

class LSRarefaction  : public DiversityCalculator  {
    
public:
    
    LSRarefaction(double c);
    
    vector<double> getValues(int ns, vector<mcmcSample>& sampling);
    
    string getTag() { return "ls"; }
    
private:
    
    double coverage;
    
};

/***********************************************************************/




#endif /* lsrarefaction_hpp */
