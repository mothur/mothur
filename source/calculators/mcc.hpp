//
//  mcc.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef mcc_hpp
#define mcc_hpp

#include <calculator.h>

/***********************************************************************/

class MCC : public ClusterMetric  {
    
public:
    MCC() : ClusterMetric("mcc") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn);
    string getCitation() { return "http://www.mothur.org/wiki/MCC"; }
    
private:
    
};

/***********************************************************************/

#endif /* mcc_hpp */
