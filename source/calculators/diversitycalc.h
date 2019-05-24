//
//  diversitycalc.h
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversitycalc_h
#define diversitycalc_h

#include "mothurout.h"
#include "sabundvector.hpp"
#include "utils.hpp"

/***********************************************************************/

class DiversityCalculator {
    
public:
    DiversityCalculator(bool rs){ m = MothurOut::getInstance();  requiresSamples = rs; }
    virtual ~DiversityCalculator(){};
    
    virtual string getTag() = 0;
    virtual bool requiresSample() { return requiresSamples; }
    
protected:
    Utils util;
    MothurOut* m;
    
    bool requiresSamples;
    
    
};
/***********************************************************************/


#endif /* diversitycalc_h */
