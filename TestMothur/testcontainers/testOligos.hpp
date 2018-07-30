//
//  testOligos.hpp
//  Mothur
//
//  Created by Sarah Westcott on 7/30/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef testOligos_hpp
#define testOligos_hpp

#include "gtest/gtest.h"
#include "oligos.h"
#include "dataset.h"

class TestOligos : public Oligos {
    
public:
    
    TestOligos();
    ~TestOligos();
    
    //MothurOut* m;
    vector<string> oligosfiles; //single, paired, indexes, comboNamesTest
};



#endif /* testOligos_hpp */
