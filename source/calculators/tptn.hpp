//
//  tptn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef tptn_hpp
#define tptn_hpp

#include "calculator.h"

/***********************************************************************/

class TPTN : public ClusterMetric  {
    
public:
    TPTN() : ClusterMetric("tptn") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/TPTN"; }
    
private:
    
};

/***********************************************************************/




#endif /* tptn_hpp */
