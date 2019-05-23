//
//  sishift.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef sishift_hpp
#define sishift_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"


/***********************************************************************/

class SIShift : public DiversityCalculator  {
    
public:
    
    SIShift();
    
    vector<double> getValues(int ns, vector<mcmcSample>& sampling);
    
    
private:
    
    
};

/***********************************************************************/




#endif /* sishift_hpp */
