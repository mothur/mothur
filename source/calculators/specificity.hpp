//
//  specificity.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef specificity_hpp
#define specificity_hpp

#include <calculator.h>

/***********************************************************************/

class Specificity : public ClusterMetric  {
    
public:
    Specificity() : ClusterMetric("spec") {};
    double getValue( long long tp,  long long tn,  long long fp,  long long fn); //ignores tp, fn
    string getCitation() { return "http://www.mothur.org/wiki/Specificity"; }
    
private:
    
};

/***********************************************************************/



#endif /* specificity_hpp */
