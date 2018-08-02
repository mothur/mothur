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


/**************************************************************************************************/

