//
//  testintegrationalignseqscommand.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/9/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testintegrationalignseqscommand_hpp
#define testintegrationalignseqscommand_hpp

#include "gtest/gtest.h"
#include "command.hpp"
#include "dataset.h"

class TestAlignSeqsIntegration  {
    
    
public:
    
    TestAlignSeqsIntegration();
    ~TestAlignSeqsIntegration();
    
    string inputDir, outputDir, setDirInputs;
    vector<string> filenames;
    CurrentFile* current;
    MothurOut* m;
    
};


#endif /* testintegrationalignseqscommand_hpp */
