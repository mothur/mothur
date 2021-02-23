//
//  picrust.hpp
//  Mothur
//
//  Created by Sarah Westcott on 11/16/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#ifndef picrust_hpp
#define picrust_hpp

#include "mothurout.h"
#include "utils.hpp"
#include "phylotree.h"
#include "sharedrabundvectors.hpp"
#include "sharedrabundfloatvectors.hpp"

/**************************************************************************************************/

class Picrust {
    
public:
    Picrust(string, string); //reference, otumap
    Picrust();
    ~Picrust();
    
    void read(string, string);
    
    void setGGOTUIDs(map<string, string>&, SharedRAbundFloatVectors*&);
    void setGGOTUIDs(map<string, string>&, SharedRAbundVectors*&);
    
        
protected:
    MothurOut* m;
    Utils util;
    
    PhyloTree* phyloTree;
    map<string, string> otuMap;
    
    void readGGOtuMap(string); //fills otuMap
    
    
};

/**************************************************************************************************/


#endif /* picrust_hpp */
