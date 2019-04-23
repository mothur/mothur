//
//  igabundance.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/3/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef igabundance_hpp
#define igabundance_hpp

#include "mothurout.h"
#include "sabundvector.hpp"


/***********************************************************************/

class IGAbundance   {
    
public:
    IGAbundance(){ m = MothurOut::getInstance(); }
    double getValues(SAbundVector* rank);
    bool requiresSample() { return true; }
    
private:
    Utils util;
    MothurOut* m;
};

/***********************************************************************/



#endif /* igabundance_hpp */
