//
//  testfastqread.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/29/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "testfastqread.h"
#include "dataset.h"

/**************************************************************************************************/
TestFastqRead::TestFastqRead() {  //setup
    m = MothurOut::getInstance();
    TestFastqDataSet data;
    reads = data.getForwardFastq();
    filenames = data.getSubsetFRFastq(100);
}
/**************************************************************************************************/
TestFastqRead::~TestFastqRead() {
    for (int i = 0; i < filenames.size(); i++) { util.mothurRemove(filenames[i]); }
     //teardown
}
/**************************************************************************************************/
//Using ATGCGTCATC & 40 39 38 37 36 35 34 33 32 31
TEST(Test_Container_FastqRead, FastqReadConstructor) {
    TestFastqRead testFastq;
    
    vector<int> scores; for (int i = 31; i < 41; i++) { scores.push_back(i); }
    Sequence seq("testSeq", "ATGCGTCATC");
    QualityScores qual("testSeq", scores);
    FastqRead read(seq, qual);
    
    EXPECT_EQ(read.getSeq(), "ATGCGTCATC");
    EXPECT_EQ(read.getScores()[0], 31);
}

TEST(Test_Container_FastqRead, FastqReadFromFileConstructor) {
    TestFastqRead testFastq;
        
    ifstream in; bool ignore; string format = "illumina1.8+";
    Utils util; util.openInputFile(testFastq.filenames[0], in);
        
    FastqRead read(in, ignore, format);
    
    EXPECT_EQ(read.getSeq(), (testFastq.reads[0]).getSeq());
    EXPECT_EQ(read.getScores()[0], testFastq.reads[0].getScores()[0]);
}
/**************************************************************************************************/
