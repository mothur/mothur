//
//  testoptimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/6/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "testoptimatrix.h"
#include "distancecommand.h"
#include "dataset.h"

/**************************************************************************************************/
TestOptiMatrix::TestOptiMatrix() {  //setup
    m = MothurOut::getInstance();
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(100); //Fasta, name, group returned
    
    string inputString = "fasta=" + filenames[0];
    m->mothurOut("/******************************************/"); m->mothurOutEndLine();
    m->mothurOut("Running command: dist.seqs(" + inputString + ")"); m->mothurOutEndLine();
    m->mothurCalling = true;
    
    Command* distCommand = new DistanceCommand(inputString);
    distCommand->execute();
    
    map<string, vector<string> > outputFilenames = distCommand->getOutputFiles();
    
    delete distCommand;
    m->mothurCalling = false;
    
    columnFile = outputFilenames["column"][0];
    m->mothurOut("/******************************************/"); m->mothurOutEndLine();
    
    inputString = "fasta=" + filenames[0] + ", output=lt";
    m->mothurOut("/******************************************/"); m->mothurOutEndLine();
    m->mothurOut("Running command: dist.seqs(" + inputString + ")"); m->mothurOutEndLine();
    m->mothurCalling = true;
    
    Command* dist2Command = new DistanceCommand(inputString);
    dist2Command->execute();
    
    outputFilenames = dist2Command->getOutputFiles();
    
    delete dist2Command;
    m->mothurCalling = false;
    
    phylipFile = outputFilenames["phylip"][0];
    m->mothurOut("/******************************************/"); m->mothurOutEndLine();
    
    blastFile = m->getTestFilePath() + "bap.100.blast";
}
/**************************************************************************************************/
TestOptiMatrix::~TestOptiMatrix() {
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); } //teardown
    m->mothurRemove(columnFile);
    m->mothurRemove(phylipFile);
}
/**************************************************************************************************
TEST_CASE("Testing OptiMatrix Class") {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", 0.03, false);
    OptiMatrix pmatrix(testOMatrix.phylipFile, "", "", 0.03, false);
    
    SECTION("Testing Column Read") {
        INFO("Using First 100 sequences of final.fasta and final.names") // Only appears on a FAIL
        
        CAPTURE(matrix.print(cout)); // Displays this variable on a FAIL
        
        CHECK(matrix.print(cout) == 112); //numdists in matrix
    }
    
    SECTION("Testing Phylip Read") {
        INFO("Using First 100 sequences of final.fasta and final.names") // Only appears on a FAIL
        
        CAPTURE(pmatrix.print(cout)); // Displays this variable on a FAIL
        
        CHECK(pmatrix.print(cout) == 112); //numdists in matrix
    }
    
    /* First few rows of matrix
     12	23	44
     10	23	32	36
     16	25	33	48
     38
     22	45	52
     
    
    SECTION("Testing isClose") {
        INFO("Sequences 0 and 1") // Only appears on a FAIL
        
        CAPTURE(matrix.isClose(0, 12));
        CHECK(matrix.isClose(0, 12) == true);
        CAPTURE(matrix.isClose(0, 44));
        CHECK(matrix.isClose(0, 44) == true);
        CAPTURE(matrix.isClose(1, 23));
        CHECK(matrix.isClose(1, 23) == true);
        CAPTURE(matrix.isClose(1, 36));
        CHECK(matrix.isClose(1, 36) == true);
    }
}*/
/**************************************************************************************************/
TEST_F(TestOptiMatrix, readColumn) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", 0.03, false);
    
    EXPECT_EQ(112,(matrix.print(cout)));
}

TEST_F(TestOptiMatrix, readPhylip) {
    TestOptiMatrix testOMatrix;
    OptiMatrix pmatrix(testOMatrix.phylipFile, "", "", 0.03, false);

    EXPECT_EQ(112,(pmatrix.print(cout)));
}

TEST_F(TestOptiMatrix, readBlast) {
    TestOptiMatrix testOMatrix;
    OptiMatrix bmatrix(testOMatrix.blastFile, "", "", 0.03, false);

    EXPECT_EQ(46032,(bmatrix.print(cout)));
}
/**************************************************************************************************/
