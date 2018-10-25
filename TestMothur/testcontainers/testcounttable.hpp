//
//  testcounttable.hpp
//  Mothur
//
//  Created by Sarah Westcott on 10/25/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testcounttable_hpp
#define testcounttable_hpp

#include "gtest/gtest.h"
#include "counttable.h"
#include "dataset.h"

class TestCountTable : public CountTable {
    
public:
    
    TestCountTable();
    ~TestCountTable();
    
    MothurOut* m;
    string namefile, groupfile, fastafile, countfile;
};


#endif /* testcounttable_hpp */
