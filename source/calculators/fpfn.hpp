//
//  fpfn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fpfn_hpp
#define fpfn_hpp

#include <calculator.h>

/***********************************************************************/

class FPFN : public ClusterMetric  {
    
public:
    FPFN() : ClusterMetric("fpfn") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn);
    string getCitation() { return "http://www.mothur.org/wiki/FPFN"; }
    
private:
    
};

/***********************************************************************/


#endif /* fpfn_hpp */
