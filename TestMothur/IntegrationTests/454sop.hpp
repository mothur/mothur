//
//  454sop.hpp
//  Mothur
//
//  Created by Sarah Westcott on 8/17/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef _54sop_hpp
#define _54sop_hpp

#include "gtest/gtest.h"
#include "command.hpp"


class Test454SOP  {
    
    
public:
    
    Test454SOP();
    ~Test454SOP();
    
    string inputDir, outputDir, setDirInputs;
    CurrentFile* current;
    MothurOut* m;
    
};


#endif /* _54sop_hpp */
