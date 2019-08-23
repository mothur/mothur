//
//  ppv.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef ppv_hpp
#define ppv_hpp

#include "calculator.h"

/***********************************************************************/

class PPV : public ClusterMetric  {
    
public:
    PPV() : ClusterMetric("ppv") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/PPV"; }
    
private:
    
};

/***********************************************************************/



#endif /* ppv_hpp */
