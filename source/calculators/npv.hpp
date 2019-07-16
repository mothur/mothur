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
    double getValue( long long tp,  long long tn,  long long fp,  long long fn);
    string getCitation() { return "http://www.mothur.org/wiki/NPV"; }
    
private:
    
};

/***********************************************************************/




#endif /* npv_hpp */
