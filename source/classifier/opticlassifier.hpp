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

/**************************************************************************************************/

class OptiClassifier : public Classify {
    
public:
    OptiClassifier() {}
    ~OptiClassifier() {}
    
    string getTaxonomy(Sequence*, string&, bool&);
    
private:
    
    
};

/**************************************************************************************************/


#endif /* opticlassifier_hpp */
