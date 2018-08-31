//
//  testOligos.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/30/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testOligos.hpp"

/**************************************************************************************************/
TestOligos::TestOligos() {  //setup
    //m = MothurOut::getInstance();
    TestDataSet data;
    oligosfiles = data.getOligosFiles(); //single, paired, indexes, comboNamesTest
}
/**************************************************************************************************/
TestOligos::~TestOligos() {
    //teardown
}
/**************************************************************************************************/

TEST(Test_Container_Oligos, Constructors) {
    TestOligos test;
    
    Oligos oligos;
    EXPECT_EQ(oligos.hasPairedPrimers(), false);
    EXPECT_EQ(oligos.hasPairedBarcodes(), false);
    
    Oligos singleOligos(test.oligosfiles[0]);
    EXPECT_EQ(singleOligos.hasPairedPrimers(), false);
    EXPECT_EQ(singleOligos.hasPairedBarcodes(), false);
    
    //did it read properly
    map<string, int> singleBarcodes = singleOligos.getBarcodes();
    int F003D000Index = singleBarcodes["AATGGTAC"];
    EXPECT_EQ("F003D000", singleOligos.getBarcodeName(F003D000Index));
    
    int MOCKGQY1XT001Index = singleBarcodes["AACCGTGTC"];
    EXPECT_EQ("MOCK.GQY1XT001", singleOligos.getBarcodeName(MOCKGQY1XT001Index));
    
    //read with reverseCompliment of reverse primer or barcode
    Oligos pairedOligos(test.oligosfiles[1]);
    EXPECT_EQ(pairedOligos.hasPairedPrimers(), true);
    EXPECT_EQ(pairedOligos.hasPairedBarcodes(), true);
    
    map<int, oligosPair> pairedBarcodes = pairedOligos.getPairedBarcodes();
    map<int, oligosPair> pairedPrimers = pairedOligos.getPairedPrimers();
    
    oligosPair F01R2A = pairedBarcodes[0];
    EXPECT_EQ("F01R2A", pairedOligos.getBarcodeName(0));
    EXPECT_EQ("CCAAC", F01R2A.forward);
    EXPECT_EQ("CAGTG", F01R2A.reverse);
    
    oligosPair V3 = pairedPrimers[0];
    EXPECT_EQ("V3", pairedOligos.getPrimerName(0));
    EXPECT_EQ("CCTACGGGAGGCAGCAG", V3.forward);
    EXPECT_EQ("CCAGCAGCCGCGGTAAT", V3.reverse);
    
    
    //read WITHOUT reverseCompliment of reverse primer or barcode
    Oligos pairedOligosNoReverse; pairedOligosNoReverse.read(test.oligosfiles[1], false);
    EXPECT_EQ(pairedOligosNoReverse.hasPairedPrimers(), true);
    EXPECT_EQ(pairedOligosNoReverse.hasPairedBarcodes(), true);
    
    pairedBarcodes = pairedOligosNoReverse.getPairedBarcodes();
    pairedPrimers = pairedOligosNoReverse.getPairedPrimers();
    
    F01R2A = pairedBarcodes[0];
    EXPECT_EQ("F01R2A", pairedOligosNoReverse.getBarcodeName(0));
    EXPECT_EQ("CCAAC", F01R2A.forward);
    EXPECT_EQ("CACTG", F01R2A.reverse);
    
    V3 = pairedPrimers[0];
    EXPECT_EQ("V3", pairedOligosNoReverse.getPrimerName(0));
    EXPECT_EQ("CCTACGGGAGGCAGCAG", V3.forward);
    EXPECT_EQ("ATTACCGCGGCTGCTGG", V3.reverse);
    
    //oligos for indexed barcode files
    Oligos indexedOligos(test.oligosfiles[2]);
    EXPECT_EQ(indexedOligos.hasPairedPrimers(), true);
    EXPECT_EQ(indexedOligos.hasPairedBarcodes(), true);
    
    pairedBarcodes = indexedOligos.getPairedBarcodes();
    pairedPrimers = indexedOligos.getPairedPrimers();
    
    oligosPair Mock3 = pairedBarcodes[0];
    EXPECT_EQ("Mock3", indexedOligos.getBarcodeName(0));
    EXPECT_EQ("NONE", Mock3.forward);
    EXPECT_EQ("CAGCTCATCAGC", Mock3.reverse);
    
    oligosPair testPrimer = pairedPrimers[0];
    EXPECT_EQ("testPrimer", indexedOligos.getPrimerName(0));
    EXPECT_EQ("NONE", testPrimer.forward);
    EXPECT_EQ("ACTYAAAKGAATTGACGG", testPrimer.reverse);
}

TEST(Test_Container_Oligos, testComboNames) {
    TestOligos test;
    
    Oligos pairedOligos(test.oligosfiles[1]);
    
    EXPECT_EQ("F01R2A.V3", pairedOligos.getGroupName(0,0));
    EXPECT_EQ("F01R2D.V5", pairedOligos.getGroupName(3,1));
    
    Oligos singleOligos(test.oligosfiles[0]);
    
    EXPECT_EQ("F003D000", singleOligos.getGroupName(0,0));
    EXPECT_EQ("F003D006", singleOligos.getGroupName(3,1));

    Oligos indexedOligos(test.oligosfiles[2]);
    EXPECT_EQ("Mock3.testPrimer", indexedOligos.getGroupName(0,0));
    EXPECT_EQ("CKD_f31.testPrimer2", indexedOligos.getGroupName(1,1));
    EXPECT_EQ("CKD_f31.testPrimer", indexedOligos.getGroupName(1,0));
    EXPECT_EQ("Mock3.testPrimer2", indexedOligos.getGroupName(0,1));
}
/**************************************************************************************************/
