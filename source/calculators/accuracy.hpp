//
//  accuracy.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef accuracy_hpp
#define accuracy_hpp

#include "calculator.h"

/***********************************************************************/

class Accuracy : public ClusterMetric  {
    
public:
    Accuracy() : ClusterMetric("accuracy") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn);
    string getCitation() { return "http://www.mothur.org/wiki/Accuracy"; }
    
private:
    
};

/***********************************************************************/


#endif /* accuracy_hpp */
