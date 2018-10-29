//
//  testsequence.h
//  Mothur
//
//  Created by Sarah Westcott on 7/27/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#ifndef Mothur_testsequence_h
#define Mothur_testsequence_h

#include "gtest/gtest.h"
#include "sequence.hpp"

class TestSequence : public Sequence {
    
    public:
    
    TestSequence();
    ~TestSequence();
    
    using Sequence::initialize;
    
    MothurOut* m;
    string fastafile;
    
    
};


#endif
