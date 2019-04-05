//
//  erarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef erarefaction_hpp
#define erarefaction_hpp

#include "mothurout.h"
#include "sabundvector.hpp"

/***********************************************************************/

class ERarefaction   {
    
public:
    ERarefaction(){ m = MothurOut::getInstance(); }
    double getValues(SAbundVector* rank, int n);
    
    
private:
    Utils util;
    MothurOut* m;
};

/***********************************************************************/





#endif /* erarefaction_hpp */
