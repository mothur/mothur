//
//  lnshift.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/14/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef lnshift_hpp
#define lnshift_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class LNShift : public DiversityCalculator  {
    
public:
    
    LNShift();
    
    vector<double> getValues(int ns, vector<mcmcSample>& sampling);
    
private:
   
    
};

/***********************************************************************/



#endif /* lnshift_hpp */
