//
//  fn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fn_hpp
#define fn_hpp

#include "calculator.h"

/***********************************************************************/

class FN : public ClusterMetric  {
    
public:
    FN() : ClusterMetric("fn") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/FN"; }
    
private:
    
};

/***********************************************************************/



#endif /* fn_hpp */
