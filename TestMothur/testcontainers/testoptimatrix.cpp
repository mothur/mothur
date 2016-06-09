//
//  testoptimatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/6/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "catch.hpp"
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
}
/**************************************************************************************************/
TestOptiMatrix::~TestOptiMatrix() {
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); } //teardown
    m->mothurRemove(columnFile);
}
/**************************************************************************************************/
TEST_CASE("Testing OptiMatrix Class") {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", 0.03, false);
    
    SECTION("Testing Read") {
        INFO("Using First 100 sequences of final.fasta and final.names") // Only appears on a FAIL
        
        CAPTURE(matrix.print(cout)); // Displays this variable on a FAIL
        
        CHECK(matrix.print(cout) == 112); //numdists in matrix
    }
    
    /* First few rows of matrix
     12	23	44
     10	23	32	36
     16	25	33	48
     38
     22	45	52
     */
    
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
}
/**************************************************************************************************/
