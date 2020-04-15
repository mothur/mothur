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
    
    ERarefaction(int inc);
    void getValues(SAbundVector* rank, vector<double>&);
    
    string getTag() { return "e"; }
    
private:
    
    int increment;
};

/***********************************************************************/





#endif /* erarefaction_hpp */
