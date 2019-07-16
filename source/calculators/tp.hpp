//
//  tp.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef tp_hpp
#define tp_hpp

#include "calculator.h"

/***********************************************************************/

class TP : public ClusterMetric  {
    
public:
    TP() : ClusterMetric("tp") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn);
    string getCitation() { return "http://www.mothur.org/wiki/TP"; }
    
private:
    
};

/***********************************************************************/



#endif /* tp_hpp */
