//
//  testfastqread.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/29/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "catch.hpp"
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
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); }
     //teardown
}
/**************************************************************************************************/

TEST_CASE("Testing FastqRead Class") {
    /*
    TestFastqRead testFastq;
    
    SECTION("Testing (Sequence, Quality) Constructor ") {
        INFO("Using ATGCGTCATC & 40 39 38 37 36 35 34 33 32 31") // Only appears on a FAIL
        
        vector<int> scores; for (int i = 31; i < 41; i++) { scores.push_back(i); }
        Sequence seq("testSeq", "ATGCGTCATC");
        QualityScores qual("testSeq", scores);
        FastqRead read(seq, qual);
        
        CAPTURE(read.getSeq()); // Displays this variable on a FAIL
        
        CHECK(read.getSeq() == "ATGCGTCATC");
        
        CAPTURE(read.getScores()[0]); // Displays this variable on a FAIL
        
        CHECK(read.getScores()[0] == 31);
    }
    
    SECTION("Testing Fastq Read From File Constructor ") {
        INFO("Using first first read in F8D0") // Only appears on a FAIL
        
        ifstream in; bool ignore; string format = "illumina1.8+";
        testFastq.m->openInputFile(testFastq.filenames[0], in);
        
        FastqRead read(in, ignore, format);
        
        CAPTURE(read.getSeq()); // Displays this variable on a FAIL
        
        CHECK(read.getSeq() == (testFastq.reads[0]).getSeq());
        
        CAPTURE(read.getScores()[0]); // Displays this variable on a FAIL
        
        CHECK(read.getScores()[0] == testFastq.reads[0].getScores()[0]);
    }
    */
}
/**************************************************************************************************/
