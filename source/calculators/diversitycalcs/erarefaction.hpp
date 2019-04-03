//
//  erarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef erarefaction_hpp
#define erarefaction_hpp

#include "calculator.h"

/***********************************************************************/

class ERarefaction : public Calculator  {
    
public:
    ERarefaction() : Calculator("erarefaction", 1, false) {};
    EstOutput getValues(SAbundVector* rank);
    EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
    string getCitation() { return "ERarefaction - Determines an emprical rarefraction curve by resampling http://www.mothur.org/wiki/erarefaction"; }
private:
    
};

/***********************************************************************/





#endif /* erarefaction_hpp */
