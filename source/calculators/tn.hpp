//
//  tn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef tn_hpp
#define tn_hpp

#include "calculator.h"

/***********************************************************************/

class TN : public ClusterMetric  {
    
public:
    TN() : ClusterMetric("tn") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/TN"; }
    
private:
    
};

/***********************************************************************/




#endif /* tn_hpp */
