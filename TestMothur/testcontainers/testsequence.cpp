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
TEST(TestSequence, SequenceConstructors) {
    
    Sequence seq("testSeq", "ATGCGTCATC");
    EXPECT_EQ(seq.getAligned(), "ATGCGTCATC");
    EXPECT_EQ(seq.getName(), "testSeq");
    EXPECT_EQ(seq.getUnaligned(), "ATGCGTCATC");
    
    Sequence seq1;
    EXPECT_EQ(seq1.getAligned(), "");
    EXPECT_EQ(seq1.getName(), "");
    EXPECT_EQ(seq1.getUnaligned(), "");
    
        
}

TEST(TestSequence, getNumBases) {
    TestSequence testSeq;
    
    
}
/**************************************************************************************************/


