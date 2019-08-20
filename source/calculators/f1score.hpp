//
//  f1score.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef f1score_hpp
#define f1score_hpp

#include "calculator.h"

/***********************************************************************/

class F1Score : public ClusterMetric  {
    
public:
    F1Score() : ClusterMetric("f1score") {};
    double getValue(double tp,  double tn,  double fp,  double fn);
    string getCitation() { return "http://www.mothur.org/wiki/F1Score"; }
    
private:
    
};

/***********************************************************************/




#endif /* f1score_hpp */
