//
//  fp.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fp_hpp
#define fp_hpp

#include <calculator.h>

/***********************************************************************/

class FP : public ClusterMetric  {
    
public:
    FP() : ClusterMetric("fp") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn);
    string getCitation() { return "http://www.mothur.org/wiki/FP"; }
    
private:
    
};

/***********************************************************************/



#endif /* fp_hpp */
