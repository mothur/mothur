//
//  testtrimoligos.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/14/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//


#include "testtrimoligos.hpp"

/**************************************************************************************************/
TestTrimOligos::TestTrimOligos() {  //setup
}
/**************************************************************************************************/
TestTrimOligos::~TestTrimOligos() {
}
/**************************************************************************************************/

TEST(Test_TrimOligos, TrimOligosConstructors) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "TGGTGAACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    int groupIndex;
    vector<int> results = noDiffSingleTrim.stripBarcode(F003D150, groupIndex);
    EXPECT_EQ(9, groupIndex);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    
    //lots of diffs allowed to force a multiMatch
    TrimOligos lotsOfDiffSingleTrim(2,8,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("TGGTGAACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    int groupIndex;
    vector<int> results = noDiffSingleTrim.stripBarcode(F003D150, groupIndex);
    EXPECT_EQ(9, groupIndex);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));

    
    //Create trimoligos classes with various constructors
    
    
    //run all public strip functions
    
    
    //run private functions
    
    
}
/**************************************************************************************************/
