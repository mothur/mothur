//
//  sirarefaction.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/23/19.
//  Copyright © 2019 Schloss Lab. All rights reserved.
//

#ifndef sirarefaction_hpp
#define sirarefaction_hpp

#include "diversityutils.hpp"
#include "diversitycalc.h"

/***********************************************************************/

class SIRarefaction : public DiversityCalculator  {
    
public:
    
    SIRarefaction(double c); // : coverage(c) { m = MothurOut::getInstance(); }
    
    vector<double> getValues(int ns, vector<mcmcSample>& sampling);
    
    string getTag() { return "si"; }
    
private:
    
    double coverage;
    
};

/***********************************************************************/




#endif /* sirarefaction_hpp */
