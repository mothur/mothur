//
//  miseqsop.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/6/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef miseqsop_hpp
#define miseqsop_hpp


#include "gtest/gtest.h"
#include "command.hpp"


class TestMiSeqSOP  {
    
    
public:
    
    TestMiSeqSOP();
    ~TestMiSeqSOP();
    
    string inputDir, outputDir, setDirInputs;
    CurrentFile* current;
    MothurOut* m;
    
};


#endif /* miseqsop_hpp */
