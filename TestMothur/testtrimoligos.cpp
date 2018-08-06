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

TEST(Test_TrimOligos, SingleDirectionStripBarcodes) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "TGGTGAACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    int groupIndex;
    vector<int> results = noDiffSingleTrim.stripBarcode(F003D150, groupIndex);
    EXPECT_EQ(9, groupIndex);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    EXPECT_EQ("CCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned()); //TGGTGAAC barcode removed
    
    //1 barcode diff, 2 primer diffs
    TrimOligos someDiffsSingleTrim(2,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("TGGTGTACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    
    results = someDiffsSingleTrim.stripBarcode(F003D150, groupIndex);
    EXPECT_EQ(9, groupIndex);
    EXPECT_EQ(1, results[0]);
    EXPECT_EQ("match", someDiffsSingleTrim.getCodeValue(results[1], 1));
    EXPECT_EQ("CCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned()); //TGGTGAAC barcode removed
    
    //nomatch option
    F003D150.setAligned("GTACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    results = someDiffsSingleTrim.stripBarcode(F003D150, groupIndex);
    EXPECT_EQ(3, results[0]);
    EXPECT_EQ("noMatch", someDiffsSingleTrim.getCodeValue(results[1], 1));
    
    //4 barcode diff, 2 primer diffs - force multiMatch
    //barcode	TTC GT G G C	F003D004
    //barcode	TTC TT G A C	F003D006
    TrimOligos lotsOfDiffsSingleTrim(2,0,4,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    Sequence F003D006("GQY1XT001CAJGU", "TTCAATACCCGTCAATTCCTTTAAGTTTCAACCTTGCGATCGTACTCCCCAGGTGGGATACTTATTGCGTTAGCTGCGGCACGCAGGGGGTCAGTCCCCGCACACCTAGTATCCATCGTTTACAGCGTGGACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCACCGTCAGTTGTCGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCCTGCCTCTCCGATACTCAAGACCTACAGTTTCAAATGCA");
    
    results = lotsOfDiffsSingleTrim.stripBarcode(F003D006, groupIndex);
    EXPECT_EQ(3, results[0]); //three diffs to find match, but matches F003D004 and F003D006
    EXPECT_EQ("multipleMatches", lotsOfDiffsSingleTrim.getCodeValue(results[1], 4));
}

TEST(Test_TrimOligos, SingleDirectionStripQualBarcodes) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "TGGTGAACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG"); //fragment
    vector<int> qualScores(53, 38);
    QualityScores F003D150Q("GQY1XT001ASWK1", qualScores); //fragment
    
    int groupIndex;
    vector<int> results = noDiffSingleTrim.stripBarcode(F003D150, F003D150Q, groupIndex);
    EXPECT_EQ(9, groupIndex);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    EXPECT_EQ("CCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG", F003D150.getUnaligned());//TGGTGAAC barcode removed
    EXPECT_EQ(45, F003D150Q.getLength()); //barcode removed
    
    
    //1 barcode diff, 2 primer diffs
    TrimOligos someDiffsSingleTrim(2,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("TGGTGCACCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG");
    F003D150Q.setScores(qualScores);
    
    results = someDiffsSingleTrim.stripBarcode(F003D150, F003D150Q, groupIndex);
    EXPECT_EQ(9, groupIndex);
    EXPECT_EQ(1, results[0]);
    EXPECT_EQ("match", someDiffsSingleTrim.getCodeValue(results[1], 1));
    EXPECT_EQ("CCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG", F003D150.getUnaligned()); //TGGTGAAC barcode removed
    EXPECT_EQ(45, F003D150Q.getLength()); //barcode removed
    
    //nomatch option
    F003D150.setAligned("TCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG");
    
    results = someDiffsSingleTrim.stripBarcode(F003D150, groupIndex);
    EXPECT_EQ(5, results[0]);
    EXPECT_EQ("noMatch", someDiffsSingleTrim.getCodeValue(results[1], 1));
    
    //4 barcode diff, 2 primer diffs - force multiMatch TTCTTGAC
    //barcode	TTC G TG G C	F003D004
    //barcode	TTC T TG A C	F003D006
    TrimOligos lotsOfDiffsSingleTrim(2,0,2,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    Sequence F003D006("GQY1XT001CAJGU", "TTCATTGCCCGTCAATTCCTTTAAGTTTCAACCTTGCGATCGTACTCCCCAGG");
    F003D150Q.setScores(qualScores);
    
    results = lotsOfDiffsSingleTrim.stripBarcode(F003D006, F003D150Q, groupIndex);
    EXPECT_EQ(2, results[0]); //three diffs to find match, but matches F003D004 and F003D006
    EXPECT_EQ("multipleMatches", lotsOfDiffsSingleTrim.getCodeValue(results[1], 2));
}

TEST(Test_TrimOligos, SingleDirectionStripPrimers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "CCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    int groupIndex;
    vector<int> results = noDiffSingleTrim.stripForward(F003D150, groupIndex);
    EXPECT_EQ(0, groupIndex);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    
    //CCGTCAATTCMTTTRAGT primer removed
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned()); //primer removed
    
    //1 barcode diff, 2 primer diffs
    TrimOligos someDiffsSingleTrim(2,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("CCTTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    
    results = someDiffsSingleTrim.stripForward(F003D150, groupIndex);
    EXPECT_EQ(0, groupIndex);
    EXPECT_EQ(1, results[0]);
    EXPECT_EQ("match", someDiffsSingleTrim.getCodeValue(results[1], 1));
    
    //TGGTGAAC barcode removed
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned()); //barcode removed
    
    //nomatch option
    F003D150.setAligned("CCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTGCCACCCAGGGGTCAATCCCCCCGGACAGCTAGCATTCATCGTTTACTGTGCGGACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT");
    
    
    results = someDiffsSingleTrim.stripForward(F003D150, groupIndex);
    EXPECT_EQ(9, results[0]);
    EXPECT_EQ("noMatch", someDiffsSingleTrim.getCodeValue(results[1], 1));
}

TEST(Test_TrimOligos, SingleDirectionStripQualPrimers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "CCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG"); //fragment
    vector<int> qualScores(45, 38);
    QualityScores F003D150Q("GQY1XT001ASWK1", qualScores); //fragment
    
    int groupIndex;
    vector<int> results = noDiffSingleTrim.stripForward(F003D150, F003D150Q, groupIndex, false);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGG", F003D150.getUnaligned());//TGGTGAAC barcode removed
    EXPECT_EQ(27, F003D150Q.getLength()); //barcode removed
    
    
    //1 barcode diff, 2 primer diffs
    TrimOligos someDiffsSingleTrim(2,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("CGGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG");
    F003D150Q.setScores(qualScores);
    
    results = someDiffsSingleTrim.stripForward(F003D150, F003D150Q, groupIndex, false);
    EXPECT_EQ(1, results[0]);
    EXPECT_EQ("match", someDiffsSingleTrim.getCodeValue(results[1], 1));
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGG", F003D150.getUnaligned()); //TGGTGAAC barcode removed
    EXPECT_EQ(27, F003D150Q.getLength()); //barcode removed
    
    //nomatch option
    F003D150.setAligned("TTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGG");
    F003D150Q.setScores(qualScores);
    
    results = someDiffsSingleTrim.stripForward(F003D150, F003D150Q, groupIndex, false);
    EXPECT_EQ(11, results[0]);
    EXPECT_EQ("noMatch", someDiffsSingleTrim.getCodeValue(results[1], 1));
}

TEST(Test_TrimOligos, SingleDirectionStripReversePrimers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");

    vector<int> results = noDiffSingleTrim.stripReverse(F003D150);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned()); //primer removed
    
    //1 barcode diff, 2 primer diffs
    TrimOligos someDiffsSingleTrim(2,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTCC");
    
    results = someDiffsSingleTrim.stripReverse(F003D150);
    EXPECT_EQ(2, results[0]);
    EXPECT_EQ("match", someDiffsSingleTrim.getCodeValue(results[1], 2));
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned()); //barcode removed
    
    //nomatch option
    F003D150.setAligned("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCCCC");
    
    results = someDiffsSingleTrim.stripReverse(F003D150);
    EXPECT_EQ(3, results[0]);
    EXPECT_EQ("noMatch", someDiffsSingleTrim.getCodeValue(results[1], 2));
}

TEST(Test_TrimOligos, SingleDirectionStripQualReversePrimers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    vector<int> qualScores(84, 38);
    QualityScores F003D150Q("GQY1XT001ASWK1", qualScores); //fragment
    
    vector<int> results = noDiffSingleTrim.stripReverse(F003D150, F003D150Q);
    EXPECT_EQ(0, results[0]);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned());
    EXPECT_EQ(67, F003D150Q.getLength()); //reverse removed
    
    //1 barcode diff, 2 primer diffs
    TrimOligos someDiffsSingleTrim(2,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, bdiffs, primers, barcodes, revPrimers
    F003D150.setAligned("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTCC");
    F003D150Q.setScores(qualScores);
    
    results = someDiffsSingleTrim.stripReverse(F003D150, F003D150Q);
    EXPECT_EQ(2, results[0]);
    EXPECT_EQ("match", someDiffsSingleTrim.getCodeValue(results[1], 2));
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTT", F003D150.getUnaligned());
    EXPECT_EQ(67, F003D150Q.getLength()); //reverse removed
    
    //nomatch option
    F003D150.setAligned("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCCCC");
    F003D150Q.setScores(qualScores);
    
    results = someDiffsSingleTrim.stripReverse(F003D150);
    EXPECT_EQ(3, results[0]);
    EXPECT_EQ("noMatch", someDiffsSingleTrim.getCodeValue(results[1], 2));
    EXPECT_EQ(84, F003D150Q.getLength()); //reverse not removed
}

TEST(Test_TrimOligos, TwoDirectionStripBarcodes) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadPaired();
    
    //no diffs allowed
    TrimOligos noDiffTwoDirectionTrim(0,0,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    Sequence F05R2F_forward("F05R2F_forward" ,"CTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    Sequence F05R2F_reverse("F05R2F_reverse" ,"GGGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    int groupIndex;
    vector<int> results = noDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, groupIndex); //cttac	gggtt
    EXPECT_EQ(8, groupIndex);
    EXPECT_EQ(0, (results[0]+results[2]));
    EXPECT_EQ("match", noDiffTwoDirectionTrim.getCodeValue(results[3], 0));
    EXPECT_EQ("ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
    
    //diffs allowed
    TrimOligos someDiffTwoDirectionTrim(0,1,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    F05R2F_forward.setAligned("CTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("CGGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    results = someDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, groupIndex); //cttac	gggtt
    EXPECT_EQ(8, groupIndex);
    EXPECT_EQ(1, (results[0]+results[2]));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[3], 1));
    EXPECT_EQ("ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
    
    F05R2F_forward.setAligned("CATACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("GGGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    results = someDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, groupIndex); //cttac	gggtt
    EXPECT_EQ(8, groupIndex);
    EXPECT_EQ(1, (results[0]+results[2]));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[3], 1));
    EXPECT_EQ("ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed

    F05R2F_forward.setAligned("GTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("CCGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    results = someDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, groupIndex); //cttac	gggtt
    EXPECT_EQ(3, (results[0]+results[2]));
    EXPECT_EQ("noMatch", someDiffTwoDirectionTrim.getCodeValue(results[3], 1));
    EXPECT_EQ("GTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
}

TEST(Test_TrimOligos, TwoDirectionStripQualBarcodes) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadPaired();
    
    //no diffs allowed
    TrimOligos noDiffTwoDirectionTrim(0,0,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    Sequence F05R2F_forward("F05R2F_forward" ,"CTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    Sequence F05R2F_reverse("F05R2F_reverse" ,"GGGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    vector<int> qualScores(82, 38);
    QualityScores F05R2F_forwardQual("F05R2F_forward", qualScores); //fragment
    QualityScores F05R2F_reverseQual("F05R2F_reverse", qualScores); //fragment
    
    int groupIndex;
    vector<int> results = noDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex); //cttac	gggtt
    EXPECT_EQ(8, groupIndex);
    EXPECT_EQ(0, (results[0]+results[2]));
    EXPECT_EQ("match", noDiffTwoDirectionTrim.getCodeValue(results[3], 0));
    EXPECT_EQ("ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
    EXPECT_EQ(77, F05R2F_forwardQual.getLength()); //barcode removed
    EXPECT_EQ(77, F05R2F_reverseQual.getLength()); //barcode removed
    
    //diffs allowed
    TrimOligos someDiffTwoDirectionTrim(0,1,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    F05R2F_forward.setAligned("CTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("CGGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    F05R2F_forwardQual.setScores(qualScores);
    F05R2F_reverseQual.setScores(qualScores);
    
    results = someDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex); //cttac	gggtt
    EXPECT_EQ(8, groupIndex);
    EXPECT_EQ(1, (results[0]+results[2]));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[3], 1));
    EXPECT_EQ("ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
    EXPECT_EQ(77, F05R2F_forwardQual.getLength()); //barcode removed
    EXPECT_EQ(77, F05R2F_reverseQual.getLength()); //barcode removed
    
    F05R2F_forward.setAligned("CATACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("GGGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    F05R2F_forwardQual.setScores(qualScores);
    F05R2F_reverseQual.setScores(qualScores);
    
    results = someDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex); //cttac	gggtt
    EXPECT_EQ(8, groupIndex);
    EXPECT_EQ(1, (results[0]+results[2]));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[3], 1));
    EXPECT_EQ("ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
    EXPECT_EQ(77, F05R2F_forwardQual.getLength()); //barcode removed
    EXPECT_EQ(77, F05R2F_reverseQual.getLength()); //barcode removed
    
    F05R2F_forward.setAligned("GTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("CCGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    F05R2F_forwardQual.setScores(qualScores);
    F05R2F_reverseQual.setScores(qualScores);
    
    results = someDiffTwoDirectionTrim.stripBarcode(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex); //cttac	gggtt
    EXPECT_EQ(3, (results[0]+results[2]));
    EXPECT_EQ("noMatch", someDiffTwoDirectionTrim.getCodeValue(results[3], 1));
    EXPECT_EQ("GTTACATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned()); //cttac barcode removed
    EXPECT_EQ("CCGTTCCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned()); //gggtt barcode removed
    EXPECT_EQ(82, F05R2F_forwardQual.getLength()); //barcode removed
    EXPECT_EQ(82, F05R2F_reverseQual.getLength()); //barcode removed
}

TEST(Test_TrimOligos, TwoDirectionStripPrimers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadPaired();
    
    //no diffs allowed
    TrimOligos noDiffTwoDirectionTrim(0,0,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    Sequence F05R2F_forward("F05R2F_forward" ,"ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    Sequence F05R2F_reverse("F05R2F_reverse" ,"CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    int groupIndex; vector<int> results = noDiffTwoDirectionTrim.stripForward(F05R2F_forward, F05R2F_reverse, groupIndex);
    EXPECT_EQ(1, groupIndex);
    EXPECT_EQ(0, (results[0]+results[2]));
    EXPECT_EQ("match", noDiffTwoDirectionTrim.getCodeValue(results[3], 0));
    EXPECT_EQ("ACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned());
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned());
    
    //allows for 2 diffs on each end - the command code decides pass/fail
    TrimOligos someDiffTwoDirectionTrim(2,0,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    F05R2F_forward.setAligned("TTTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("GCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    results = someDiffTwoDirectionTrim.stripForward(F05R2F_forward, F05R2F_reverse, groupIndex);
    EXPECT_EQ(1, groupIndex);
    EXPECT_EQ(2, (results[0]+results[2]));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[3], 2));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[1], 2));
    EXPECT_EQ("ACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned());
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned());
    
    //allows for 2 diffs on each end - the command code decides pass/fail
    F05R2F_forward.setAligned("TTTTGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("GGCGTGAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    
    results = someDiffTwoDirectionTrim.stripForward(F05R2F_forward, F05R2F_reverse, groupIndex);
    EXPECT_EQ(1, groupIndex);
    EXPECT_EQ(5, (results[0]+results[2]));
    EXPECT_EQ("noMatch", someDiffTwoDirectionTrim.getCodeValue(results[3], 2));
    EXPECT_EQ("noMatch", someDiffTwoDirectionTrim.getCodeValue(results[1], 2));
    EXPECT_EQ("TTTTGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned());
    EXPECT_EQ("GGCGTGAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned());
}

TEST(Test_TrimOligos, TwoDirectionStripQualPrimers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadPaired();
    
    //no diffs allowed
    TrimOligos noDiffTwoDirectionTrim(0,0,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    Sequence F05R2F_forward("F05R2F_forward" ,"ATTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    vector<int> qualScores(77, 38);
    QualityScores F05R2F_forwardQual("F05R2F_forward", qualScores); //fragment
    Sequence F05R2F_reverse("F05R2F_reverse" ,"CCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    QualityScores F05R2F_reverseQual("F05R2F_reverse", qualScores); //fragment
    
    int groupIndex; vector<int> results = noDiffTwoDirectionTrim.stripForward(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex);
    EXPECT_EQ(1, groupIndex);
    EXPECT_EQ(0, (results[0]+results[2]));
    EXPECT_EQ("match", noDiffTwoDirectionTrim.getCodeValue(results[3], 0));
    EXPECT_EQ("ACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned());
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned());
    EXPECT_EQ(57, F05R2F_forwardQual.getLength()); //primer removed
    EXPECT_EQ(58, F05R2F_reverseQual.getLength()); //primer removed
    
    //allows for 2 diffs on each end - the command code decides pass/fail
    TrimOligos someDiffTwoDirectionTrim(2,0,0,0,testTrim.oligos.ipprimers, testTrim.oligos.ipbarcodes, false); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, hasIndex
    F05R2F_forward.setAligned("TTTAGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("GCCGTCAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    F05R2F_forwardQual.setScores(qualScores);
    F05R2F_reverseQual.setScores(qualScores);
    
    results = someDiffTwoDirectionTrim.stripForward(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex);
    EXPECT_EQ(1, groupIndex);
    EXPECT_EQ(2, (results[0]+results[2]));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[3], 2));
    EXPECT_EQ("match", someDiffTwoDirectionTrim.getCodeValue(results[1], 2));
    EXPECT_EQ("ACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned());
    EXPECT_EQ("TTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned());
    EXPECT_EQ(57, F05R2F_forwardQual.getLength()); //primer removed
    EXPECT_EQ(58, F05R2F_reverseQual.getLength()); //primer removed
    
    //allows for 2 diffs on each end - the command code decides pass/fail
    F05R2F_forward.setAligned("TTTTGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG");
    F05R2F_reverse.setAligned("GGCGTGAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC");
    F05R2F_forwardQual.setScores(qualScores);
    F05R2F_reverseQual.setScores(qualScores);
    
    results = someDiffTwoDirectionTrim.stripForward(F05R2F_forward, F05R2F_reverse, F05R2F_forwardQual, F05R2F_reverseQual, groupIndex);
    EXPECT_EQ(1, groupIndex);
    EXPECT_EQ(5, (results[0]+results[2]));
    EXPECT_EQ("noMatch", someDiffTwoDirectionTrim.getCodeValue(results[3], 2));
    EXPECT_EQ("noMatch", someDiffTwoDirectionTrim.getCodeValue(results[1], 2));
    EXPECT_EQ("TTTTGATACCCGGGTAGTCCACGCAGTAAACGATGCATGCTAACTGTCAGGTGCGTTGAGCGCGGTGCGATGCAGCG", F05R2F_forward.getUnaligned());
    EXPECT_EQ("GGCGTGAATTCATTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATGCTTAACGCTTTCGCTGTACCGCCC", F05R2F_reverse.getUnaligned());
    EXPECT_EQ(77, F05R2F_forwardQual.getLength()); //primer removed
    EXPECT_EQ(77, F05R2F_reverseQual.getLength()); //primer removed
}

TEST(Test_TrimOligos, SingleDirectionStripLinkers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer, testTrim.oligos.linker, testTrim.oligos.spacer); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimers, linker, spacer
    
    Sequence F003D150("GQY1XT001ASWK1", "TGACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    int result = noDiffSingleTrim.stripLinker(F003D150);
    EXPECT_EQ(0, result);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    
    //1 linker diff
    TrimOligos someDiffsSingleTrim(0,0,1,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer, testTrim.oligos.linker, testTrim.oligos.spacer);
    F003D150.setAligned("AGACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    result = someDiffsSingleTrim.stripLinker(F003D150);
    EXPECT_EQ(1, result);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    
    //nomatch option
    F003D150.setAligned("TTCATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    result = someDiffsSingleTrim.stripLinker(F003D150);
    EXPECT_EQ(7, result);
    EXPECT_EQ("TTCATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
   
}

TEST(Test_TrimOligos, SingleDirectionStripQualLinkers) {
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer, testTrim.oligos.linker, testTrim.oligos.spacer); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimers, linker, spacer
    
    Sequence F003D150("GQY1XT001ASWK1", "TGACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    vector<int> qualScores(88, 38);
    QualityScores F003D150Q("F003D150Q", qualScores); //fragment
    
    int result = noDiffSingleTrim.stripLinker(F003D150, F003D150Q);
    EXPECT_EQ(0, result);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    EXPECT_EQ(84, F003D150Q.getLength());
    
    //1 linker diff
    TrimOligos someDiffsSingleTrim(0,0,1,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer, testTrim.oligos.linker, testTrim.oligos.spacer);
    F003D150.setAligned("AGACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    F003D150Q.setScores(qualScores);
    
    result = someDiffsSingleTrim.stripLinker(F003D150, F003D150Q);
    EXPECT_EQ(1, result);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    EXPECT_EQ(84, F003D150Q.getLength());
    
    //nomatch option
    F003D150.setAligned("TTCATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    F003D150Q.setScores(qualScores);
    
    result = someDiffsSingleTrim.stripLinker(F003D150, F003D150Q);
    EXPECT_EQ(7, result);
    EXPECT_EQ("TTCATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    EXPECT_EQ(88, F003D150Q.getLength());
    
}

TEST(Test_TrimOligos, SingleDirectionStripSpacers) { //CACTG, CCAAC
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed
    TrimOligos noDiffSingleTrim(0,0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer, testTrim.oligos.linker, testTrim.oligos.spacer); //pdiffs, bdiffs, ldiffs, sdiffs, primers, barcodes, revPrimers, linker, spacer
    
    Sequence F003D150("GQY1XT001ASWK1", "CCAACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    int result = noDiffSingleTrim.stripSpacer(F003D150);
    EXPECT_EQ(0, result);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    
    //1 linker diff
    TrimOligos someDiffsSingleTrim(0,0,0,1,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer, testTrim.oligos.linker, testTrim.oligos.spacer);
    F003D150.setAligned("CAGTGTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    result = someDiffsSingleTrim.stripSpacer(F003D150);
    EXPECT_EQ(1, result);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    
    //nomatch option
    F003D150.setAligned("CGTACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    result = someDiffsSingleTrim.stripSpacer(F003D150);
    EXPECT_EQ(2, result);
    EXPECT_EQ("CGTACTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", F003D150.getUnaligned());
    
}

TEST(Test_TrimOligos, SingleDirectionFindForward) { //CCGTCAATTCMTTTRAGT
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed - unaligned
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, rpdiffs, bdiffs, primers, barcodes, revPrimers
    
    Sequence F003D150("GQY1XT001ASWK1", "CCAACCCGTCAATTCMTTTRAGTTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG");
    
    int primerStart, primerEnd;
    vector<int> results = noDiffSingleTrim.findForward(F003D150, primerStart, primerEnd);
    EXPECT_EQ(5, primerStart);
    EXPECT_EQ(23, primerEnd);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    string primerRemoved = (F003D150.getUnaligned().substr(primerEnd));
    string primerKept = (F003D150.getUnaligned().substr(primerStart));
    EXPECT_EQ("CCGTCAATTCMTTTRAGTTATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", primerKept);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", primerRemoved);
    
    //no diffs allowed - aligned gap between primer and sequence
    F003D150.setAligned("......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    Sequence F003D150_1("GQY1XT001ASWK1", "......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    Sequence F003D150_2("GQY1XT001ASWK1", "......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");

    int countBases = 0;
    string temp = F003D150.getAligned();
    map<int, int> mapAligned;
    for (int i = 0; i < temp.length(); i++) {
        if (isalpha(temp[i])) { mapAligned[countBases] = i; countBases++; } //maps location in unaligned -> location in aligned.
    }                                                   //ie. the 3rd base may be at spot 10 in the alignment
    
    results = noDiffSingleTrim.findForward(F003D150, primerStart, primerEnd);
    EXPECT_EQ(5, primerStart);
    EXPECT_EQ(23, primerEnd);
    F003D150_1.filterToPos(mapAligned[primerEnd-1]+1);
    string primerRemovedKeepDots = F003D150_1.getAligned();
    primerRemoved = F003D150.getAligned().substr(mapAligned[primerEnd-1]+1);
    primerKept = F003D150.getAligned().substr(mapAligned[primerStart]);
    F003D150_2.filterToPos(mapAligned[primerStart]);
    string primerKeptKeepDots = F003D150_2.getAligned();
    
    EXPECT_EQ("CCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerKept);
    EXPECT_EQ("-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerRemoved);
    EXPECT_EQ(".............CCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerKeptKeepDots);
    EXPECT_EQ("................................TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerRemovedKeepDots);
    
    //no diffs allowed - aligned no gap between primer and sequence
    F003D150.setAligned("......CC--AACCCGTCAATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    F003D150_1.setAligned("......CC--AACCCGTCAATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    F003D150_2.setAligned("......CC--AACCCGTCAATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    
    countBases = 0; temp = F003D150.getAligned(); mapAligned.clear();
    for (int i = 0; i < temp.length(); i++) { if (isalpha(temp[i])) { mapAligned[countBases] = i; countBases++; }  }
    
    results = noDiffSingleTrim.findForward(F003D150, primerStart, primerEnd);
    EXPECT_EQ(5, primerStart);
    EXPECT_EQ(23, primerEnd);
    F003D150_1.filterToPos(mapAligned[primerEnd-1]+1);
    primerRemovedKeepDots = F003D150_1.getAligned();
    primerRemoved = F003D150.getAligned().substr(mapAligned[primerEnd-1]+1);
    primerKept = F003D150.getAligned().substr(mapAligned[primerStart]);
    F003D150_2.filterToPos(mapAligned[primerStart]);
    primerKeptKeepDots = F003D150_2.getAligned();
    
    EXPECT_EQ("CCGTCAATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerKept);
    EXPECT_EQ("TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerRemoved);
    EXPECT_EQ(".............CCGTCAATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerKeptKeepDots);
    EXPECT_EQ("...............................TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerRemovedKeepDots);

    //2 diffs allowed //CCGTCAATTCMTTTRAGT
    TrimOligos someDiffSingleTrim(2,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, nullVector); //pdiffs, rpdiffs, bdiffs, primers, barcodes, revPrimers
                                    //CCGTCAATTCMTTTRAGT
    F003D150.setAligned("......CC--AACCCGTGTATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    F003D150_1.setAligned("......CC--AACCCGTGTATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    F003D150_2.setAligned("......CC--AACCCGTGTATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    
    results = someDiffSingleTrim.findForward(F003D150, primerStart, primerEnd);
    EXPECT_EQ(5, primerStart);
    EXPECT_EQ(23, primerEnd);
    F003D150_1.filterToPos(mapAligned[primerEnd-1]+1);
    primerRemovedKeepDots = F003D150_1.getAligned();
    primerRemoved = F003D150.getAligned().substr(mapAligned[primerEnd-1]+1);
    primerKept = F003D150.getAligned().substr(mapAligned[primerStart]);
    F003D150_2.filterToPos(mapAligned[primerStart]);
    primerKeptKeepDots = F003D150_2.getAligned();
    EXPECT_EQ("CCGTGTATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerKept);
    EXPECT_EQ("TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerRemoved);
    EXPECT_EQ(".............CCGTGTATTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerKeptKeepDots);
    EXPECT_EQ("...............................TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", primerRemovedKeepDots);
                                    //CCGTCAATTCMTTTRAGT
    F003D150.setAligned("......CC--AACCCGTGTTTTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    
    results = someDiffSingleTrim.findForward(F003D150, primerStart, primerEnd);
    EXPECT_EQ(0, primerStart); //indicates failure
    EXPECT_EQ(0, primerEnd); //indicates failure
    EXPECT_EQ("......CC--AACCCGTGTTTTCMTTTRAGTTATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC--GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....", F003D150.getAligned());
}

TEST(Test_TrimOligos, SingleDirectionFindReverse) { //ATTACCGCGGCTGCTGG
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed - unaligned
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, rpdiffs, bdiffs, primers, barcodes, revPrimers
    
                                                                                                                           //ATTACCGCGGCTGCTGG
    Sequence F003D150("GQY1XT001ASWK1", "TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGGATTACCGCGGCTGCTGGTGGTCCAGT");
    
    int primerStart, primerEnd;
    vector<int> results = noDiffSingleTrim.findReverse(F003D150, primerStart, primerEnd);
    EXPECT_EQ(84, primerStart);
    EXPECT_EQ(101, primerEnd);
    EXPECT_EQ("match", noDiffSingleTrim.getCodeValue(results[1], 0));
    string primerRemoved = (F003D150.getUnaligned().substr(0, primerStart));
    string primerKept = (F003D150.getUnaligned().substr(0, primerEnd));
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGGATTACCGCGGCTGCTGG", primerKept);
    EXPECT_EQ("TATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTATTACCGCGGCTGCTGG", primerRemoved);
    
    //no diffs allowed - aligned gap between primer and sequence                                                      ATTACCGCGGCTGCTGG
    F003D150.setAligned("......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACCGCGGCTGCTGG-GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    Sequence F003D150_1("GQY1XT001ASWK1", "......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACCGCGGCTGCTGG-GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    Sequence F003D150_2("GQY1XT001ASWK1", "......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACCGCGGCTGCTGG-GAGCCCGGCAGTTTATTACC-GCGGC-GCTGG.....");
    
    int countBases = 0;
    string temp = F003D150.getAligned();
    map<int, int> mapAligned;
    for (int i = 0; i < temp.length(); i++) {
        if (isalpha(temp[i])) { mapAligned[countBases] = i; countBases++; } //maps location in unaligned -> location in aligned.
    }                                                   //ie. the 3rd base may be at spot 10 in the alignment
    
    results = noDiffSingleTrim.findReverse(F003D150, primerStart, primerEnd);
    EXPECT_EQ(76, primerStart);
    EXPECT_EQ(93, primerEnd);
    F003D150_1.filterFromPos(mapAligned[primerStart]);
    string primerRemovedKeepDots = F003D150_1.getAligned();
    primerRemoved = F003D150.getAligned().substr(0, mapAligned[primerStart]);
    primerKept = F003D150.getAligned().substr(0, mapAligned[primerEnd-1]+1);
    F003D150_2.filterFromPos(mapAligned[primerEnd-1]+1);
    string primerKeptKeepDots = F003D150_2.getAligned();
    
    EXPECT_EQ("......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACCGCGGCTGCTGG", primerKept);
    EXPECT_EQ("......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACCGCGGCTGCTGG......................................", primerKeptKeepDots);
    EXPECT_EQ("......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-", primerRemoved);
    EXPECT_EQ("......CC--AACCCGTCAATTCMTTTRAGT-TATCTATGC--ATTTCACCG-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-.......................................................", primerRemovedKeepDots);
     
    //2 diffs allowed
    TrimOligos someDiffSingleTrim(0,2,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, bdiffs, primers, barcodes, revPrimers
                                                                                                                    //ATTACCGCGGCTGCTGG
    F003D150.setAligned("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTGG-TTACG.....");
    F003D150_1.setAligned("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTGG-TTACG.....");
    F003D150_2.setAligned("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTGG-TTACG.....");
    
    mapAligned.clear(); countBases=0;
    temp = F003D150.getAligned();
    for (int i = 0; i < temp.length(); i++) {
        if (isalpha(temp[i])) { mapAligned[countBases] = i; countBases++; } //maps location in unaligned -> location in aligned.
    }
    
    results = someDiffSingleTrim.findReverse(F003D150, primerStart, primerEnd);
    EXPECT_EQ(75, primerStart);
    EXPECT_EQ(92, primerEnd);
    F003D150_1.filterFromPos(mapAligned[primerStart]);
    primerRemovedKeepDots = F003D150_1.getAligned();
    primerRemoved = F003D150.getAligned().substr(0, mapAligned[primerStart]);
    primerKept = F003D150.getAligned().substr(0, mapAligned[primerEnd-1]+1);
    F003D150_2.filterFromPos(mapAligned[primerEnd-1]+1);
    primerKeptKeepDots = F003D150_2.getAligned();
    
    EXPECT_EQ("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTGG", primerKept);
    EXPECT_EQ("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTGG...........", primerKeptKeepDots);
    EXPECT_EQ("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-", primerRemoved);
    EXPECT_EQ("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-............................", primerRemovedKeepDots);
    
    //no match option
    F003D150.setAligned("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTAA-TTACG.....");
    F003D150_1.setAligned("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTAA-TTACG.....");
    F003D150_2.setAligned("......CC--AA-C-T-ACACCACGCATTCC--GCATACTTCTCGCCCACTC-ATTACC-TTGCS-GCGGCTGCTGG-GAGCCCGGCAGTTT-ATTAGGGCGGCTGCTAA-TTACG.....");
    
    mapAligned.clear(); countBases=0;
    temp = F003D150.getAligned();
    for (int i = 0; i < temp.length(); i++) {
        if (isalpha(temp[i])) { mapAligned[countBases] = i; countBases++; } //maps location in unaligned -> location in aligned.
    }
    
    results = someDiffSingleTrim.findReverse(F003D150, primerStart, primerEnd);
    EXPECT_EQ(0, primerStart);
    EXPECT_EQ(0, primerEnd);
}

TEST(Test_TrimOligos, SingleDirectionReverseOligos) { //ATTACCGCGGCTGCTGG
    TestTrimOligos testTrim;
    testTrim.oligos.loadSingle();
    
    //no diffs allowed - unaligned
    TrimOligos noDiffSingleTrim(0,0,0,testTrim.oligos.primers, testTrim.oligos.barcodes, testTrim.oligos.revPrimer); //pdiffs, rpdiffs, bdiffs, primers, barcodes, revPrimers


    string testOligos = "ATTACCGCGGCTGCTGG";
    EXPECT_EQ("CCAGCAGCCGCGGTAAT"  , noDiffSingleTrim.reverseOligo(testOligos));
    
    testOligos = "CCGTCAATTCMTTTRAGT";
    EXPECT_EQ("ACTYAAAKGAATTGACGG"  , noDiffSingleTrim.reverseOligo(testOligos));
    
}


/**************************************************************************************************/

