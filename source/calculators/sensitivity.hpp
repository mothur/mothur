//
//  sensitivity.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef sensitivity_hpp
#define sensitivity_hpp

#include "calculator.h"

/***********************************************************************/

class Sensitivity : public ClusterMetric  {
    
public:
    Sensitivity() : ClusterMetric("sens") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn); //ignores tn, fp
    string getCitation() { return "http://www.mothur.org/wiki/Sensitivity"; }
    
private:
    
};

/***********************************************************************/

#endif /* sensitivity_hpp */
