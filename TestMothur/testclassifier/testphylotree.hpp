//
//  testphylotree.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/29/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//

#ifndef testphylotree_hpp
#define testphylotree_hpp

#include "phylotree.h"

class TestPhyloTree : public PhyloTree {
    
    
public:
    
    TestPhyloTree();
    ~TestPhyloTree();
    
    MothurOut* m;
    
    PhyloTree phylo;
    
    //using PhyloTree::
    //using PhyloTree::  
    
};


#endif /* testphylotree_hpp */
