//
//  erarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef erarefaction_hpp
#define erarefaction_hpp

#include "diversitycalc.h"

/***********************************************************************/

class ERarefaction : public DiversityCalculator   {
    
public:
    
    ERarefaction();
    double getValues(SAbundVector* rank, int n);
    
    
    
private:
    
};

/***********************************************************************/





#endif /* erarefaction_hpp */
