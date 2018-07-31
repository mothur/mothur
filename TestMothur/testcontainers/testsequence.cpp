//
//  testsequence.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/23/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//


#include "testsequence.h"
#include "dataset.h"

/**************************************************************************************************/
TestSequence::TestSequence() {  //setup
    m = MothurOut::getInstance();
    
    TestDataSet data;
    vector<string> filenames = data.getSubsetFNGFiles();
    fastafile = filenames[0];
    
}
/**************************************************************************************************/
TestSequence::~TestSequence() {}//teardown
/**************************************************************************************************/
TEST(TestContainer_Sequence, SequenceConstructors) {
    
    Sequence seq("testSeq", "ATGCGTCATC");
    EXPECT_EQ(seq.getAligned(), "ATGCGTCATC");
    EXPECT_EQ(seq.getName(), "testSeq");
    EXPECT_EQ(seq.getUnaligned(), "ATGCGTCATC");
    
    Sequence seq1;
    EXPECT_EQ(seq1.getAligned(), "");
    EXPECT_EQ(seq1.getName(), "");
    EXPECT_EQ(seq1.getUnaligned(), "");
    
        
}

TEST(TestContainer_Sequence, setGets) {
    
    Sequence seq; seq.setName("mothurSeq");
    EXPECT_EQ("mothurSeq", seq.getName());
    
    seq.setAligned("A--TGC-G-TCA--TC");
    EXPECT_EQ("A--TGC-G-TCA--TC", seq.getAligned());
    
    seq.setUnaligned("A--TGC-G-TCA--TC");
    EXPECT_EQ("ATGCGTCATC", seq.getUnaligned());
    
    seq.reverseComplement();
    EXPECT_EQ("GATGACGCAT", seq.getUnaligned());
    
    seq.setComment("This is my sequence comment. It can contain anything I want including numbers and symbols & % -- 234");
    EXPECT_EQ("This is my sequence comment. It can contain anything I want including numbers and symbols & % -- 234", seq.getComment());
    
    Sequence seq2("testSeq", "ATNNGTCATC");
    EXPECT_EQ("testSeq\tATNNGTCATC", seq2.getInlineSeq());
    
    EXPECT_EQ(2, seq2.getNumNs());
    EXPECT_EQ(0, seq.getNumNs());
    
    Sequence seq3("testSeq", "ATGCGTCATC");
    EXPECT_EQ(10, seq3.getNumBases());
    EXPECT_EQ(1, seq3.getStartPos());
    EXPECT_EQ(10, seq3.getEndPos());
    
    Sequence seq4("testSeq", "..A--TGC-G-TCA--TC..");
    EXPECT_EQ(3, seq4.getStartPos());
    EXPECT_EQ(18, seq4.getEndPos());
    
    seq3.trim(5);
    EXPECT_EQ("ATGCG", seq3.getAligned());
    
    seq4.trim(5);
    EXPECT_EQ("ATGCG", seq4.getAligned());
    
    seq4.padToPos(3);
    EXPECT_EQ("..GCG", seq4.getAligned());
    
    seq4.padFromPos(3);
    EXPECT_EQ("..G..", seq4.getAligned());
    
    seq.setAligned("ATGCG");
    seq.filterToPos(2);
    EXPECT_EQ("..GCG", seq.getAligned());
    
    seq.setAligned("ATGCG");
    seq.filterFromPos(2);
    EXPECT_EQ("AT...", seq.getAligned());
    
    seq.setAligned("..GCG");
    seq.filterFromPos(2);
    EXPECT_EQ(".....", seq.getAligned());
    
    seq.setAligned("..A--TGC-G-TCA--TC..");
    EXPECT_EQ(20, seq.getAlignLength());
    EXPECT_EQ(0, seq.getAmbigBases());
    
    seq.setAligned("..A--MGC-G-TCA--TC..");
    EXPECT_EQ(1, seq.getAmbigBases());
    
    seq.removeAmbigBases();
    EXPECT_EQ(0, seq.getAmbigBases());
    EXPECT_EQ("..A---GC-G-TCA--TC..", seq.getAligned());
    
    EXPECT_EQ(1, seq.getLongHomoPolymer());
    seq.setAligned("..A--AAA-G-TCA--TC..");
    EXPECT_EQ(4, seq.getLongHomoPolymer());
    
    EXPECT_EQ("0000231031", seq.convert2ints());
}

/**************************************************************************************************/


