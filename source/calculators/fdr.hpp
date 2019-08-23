//
//  fdr.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fdr_hpp
#define fdr_hpp

#include "calculator.h"

/***********************************************************************/

class FDR : public ClusterMetric  {
    
public:
    FDR() : ClusterMetric("fdr") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/FDR"; }
    
private:
    
};

/***********************************************************************/



#endif /* fdr_hpp */
