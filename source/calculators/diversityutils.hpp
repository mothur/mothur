//
//  diversityutils.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef diversityutils_hpp
#define diversityutils_hpp

#include "mothurout.h"
#include "utils.hpp"

/***********************************************************************/

class DiversityUtils   {
    
public:
    DiversityUtils(){ m = MothurOut::getInstance(); }
    
    
    
private:
    Utils util;
    MothurOut* m;
};

/***********************************************************************/



#endif /* diversityutils_hpp */
