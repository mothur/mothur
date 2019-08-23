//
//  fpfn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fpfn_hpp
#define fpfn_hpp

#include "calculator.h"

/***********************************************************************/

class FPFN : public ClusterMetric  {
    
public:
    FPFN() : ClusterMetric("fpfn") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/FPFN"; }
    
private:
    
};

/***********************************************************************/


#endif /* fpfn_hpp */
