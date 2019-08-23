//
//  npv.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef npv_hpp
#define npv_hpp

#include "calculator.h"

/***********************************************************************/

class NPV : public ClusterMetric  {
    
public:
    NPV() : ClusterMetric("npv") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/NPV"; }
    
private:
    
};

/***********************************************************************/




#endif /* npv_hpp */
