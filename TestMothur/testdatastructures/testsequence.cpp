//
//  testsequence.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/23/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "catch.hpp"
#include "sequence.hpp"


TEST_CASE("Testing Sequence Class") {
    Sequence seq;
    
    SECTION("test constructor - string, string") {
        INFO("Using TestSeq, atgcatgc") // Only appears on a FAIL
        Sequence seq2("TestSeq", "atgcatgc");
        CAPTURE(seq2.getInlineSeq()); // Displays this variable on a FAIL
        
        CHECK(seq2.getInlineSeq() == "TestSeq\tatgcatgc");
    }
    
    SECTION("setting / getting name") {
        INFO("Using TestSeq") // Only appears on a FAIL
        seq.setName("TestSeq");
        CAPTURE(seq.getName()); // Displays this variable on a FAIL
        
        CHECK(seq.getName() == "TestSeq");
    }
    
    SECTION("test setAligned / get Aligned") {
        INFO("Using ..atgc--atgc..") // Only appears on a FAIL
        seq.setAligned("..atgc--atgc..");
        CAPTURE(seq.getAligned()); // Displays this variable on a FAIL
        
        CHECK(seq.getAligned() == "..atgc--atgc..");
    }
    
    SECTION("test setUnaligned / getUnaligned") {
        INFO("Using ..atgc--atgc..") // Only appears on a FAIL
        seq.setUnaligned("..atgc--atgc..");
        CAPTURE(seq.getUnaligned()); // Displays this variable on a FAIL
        
        CHECK(seq.getUnaligned() == "atgcatgc");
    }
    
   //more tests need to be added - just a start to set up testing project and model
}

