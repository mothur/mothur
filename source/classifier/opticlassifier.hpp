//
//  opticlassifier.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/12/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef opticlassifier_hpp
#define opticlassifier_hpp

#include "mothur.h"
#include "classify.h"
#include "optidata.hpp"

/**************************************************************************************************/

class OptiClassifier : public Classify {
    
public:
    OptiClassifier(string reffasta, string reftax, string mothurVersion);
    ~OptiClassifier() {}
    
    string getTaxonomy(Sequence*, string&, bool&) { return "not done yet"; }
    
private:
    OptiData* matrix;
    
        
    vector< vector<string> > binReferences();

    
};

/**************************************************************************************************/


#endif /* opticlassifier_hpp */
