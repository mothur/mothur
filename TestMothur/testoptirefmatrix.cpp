//
//  testoptirefmatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/24/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testoptirefmatrix.hpp"
#include "dataset.h"

/**************************************************************************************************/
TestOptiRefMatrix::TestOptiRefMatrix() {  //setup
    m = MothurOut::getInstance();
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(); //Fasta, name, group returned
    columnFile = data.getSubsetFNGDistFile();
    phylipFile = data.getSubsetFNGPhylipDistFile();
    reffilenames = data.getOptiRefFiles(); //fasta, count, column, phylip, list, betweendist returned
}
/**************************************************************************************************/
TestOptiRefMatrix::~TestOptiRefMatrix() {}
/**************************************************************************************************/

//distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
TEST(TestOptiRefMatrix, readColumnDenovo) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
    
    //EXPECT_EQ(160,(matrix.print(cout)));
    EXPECT_EQ(160,(matrix.getNumDists()));
    EXPECT_EQ(34,(matrix.getNumFitDists()));
    EXPECT_EQ(54,(matrix.getNumRefDists()));
    
    vector<long long> refSeqs = matrix.getRefSeqs();
    string Expected_ReturnResults = "1357111516181920222527283132333435373840414445474851535556";
    string ReturnResults = "";
    for (long long i = 0; i < refSeqs.size(); i++) { ReturnResults += toString(refSeqs[i]);  }

    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    long long sanityCheck = matrix.getNumDists() - (matrix.getNumFitDists() + matrix.getNumRefDists());
    EXPECT_EQ(72,sanityCheck); //number of inbetween dists
}

//distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
TEST(TestOptiRefMatrix, readPhylipDenovo) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.phylipFile, "phylip", testOMatrix.filenames[1], "name", 0.03, 50);
    
    //EXPECT_EQ(160,(matrix.print(cout)));
    EXPECT_EQ(160,(matrix.getNumDists()));
    EXPECT_EQ(34,(matrix.getNumFitDists()));
    EXPECT_EQ(44,(matrix.getNumRefDists()));
    
    vector<long long> refSeqs = matrix.getRefSeqs();
    string Expected_ReturnResults = "136791014161819212224282930323335363839414344455152535556";
    string ReturnResults = "";
    for (long long i = 0; i < refSeqs.size(); i++) { ReturnResults += toString(refSeqs[i]);  }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    long long sanityCheck = matrix.getNumDists() - (matrix.getNumFitDists() + matrix.getNumRefDists());
    EXPECT_EQ(82,sanityCheck); //number of inbetween dists
}

//refdistfile, refname or refcount, refformat, refdistformat, cutoff, fitdistfile, fitname or fitcount, fitformat, fitdistformat, betweendistfile, betweendistformat - files for reference
TEST(TestOptiRefMatrix, readColumnReference) {
    TestOptiRefMatrix testOMatrix;
    
    OptiRefMatrix matrix(testOMatrix.reffilenames[2], testOMatrix.reffilenames[1], "count", "column", 0.03, testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", testOMatrix.reffilenames[5], "column");
    
    //EXPECT_EQ(113772,(matrix.print(cout)));
    EXPECT_EQ(56886,(matrix.getNumDists())); //unique dists 56886*2=113772
    EXPECT_EQ(80,(matrix.getNumFitDists()));
    EXPECT_EQ(56675,(matrix.getNumRefDists()));
    
    long long sanityCheck = matrix.getNumDists() - (matrix.getNumFitDists() + matrix.getNumRefDists());
    EXPECT_EQ(131,sanityCheck); //number of inbetween dists
}

//refdistfile, refname or refcount, refformat, refdistformat, cutoff, fitdistfile, fitname or fitcount, fitformat, fitdistformat, betweendistfile, betweendistformat - files for reference
TEST(TestOptiRefMatrix, readPhylipReference) {
    TestOptiRefMatrix testOMatrix;
    
    OptiRefMatrix matrix(testOMatrix.reffilenames[3], testOMatrix.reffilenames[1], "count", "phylip", 0.03, testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", testOMatrix.reffilenames[5], "column");
    
    //EXPECT_EQ(113772,(matrix.print(cout)));
    EXPECT_EQ(56886,(matrix.getNumDists())); //unique dists 56886*2=113772
    EXPECT_EQ(80,(matrix.getNumFitDists()));
    EXPECT_EQ(56675,(matrix.getNumRefDists()));
    
    long long sanityCheck = matrix.getNumDists() - (matrix.getNumFitDists() + matrix.getNumRefDists());
    EXPECT_EQ(131,sanityCheck); //number of inbetween dists
}

TEST(TestOptiRefMatrix, getNumCLose) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
    
    EXPECT_EQ(1,(matrix.getNumClose(0)));
    EXPECT_EQ(2,(matrix.getNumClose(5)));
    EXPECT_EQ(3,(matrix.getNumClose(10)));
    EXPECT_EQ(7,(matrix.getNumClose(15)));
    EXPECT_EQ(2,(matrix.getNumClose(20)));
    
    EXPECT_EQ(1,(matrix.getNumFitClose(0)));
    EXPECT_EQ(0,(matrix.getNumFitClose(5)));
    EXPECT_EQ(1,(matrix.getNumFitClose(10)));
    EXPECT_EQ(4,(matrix.getNumFitClose(15)));
    EXPECT_EQ(0,(matrix.getNumFitClose(20)));
    
    EXPECT_EQ(0,(matrix.getNumRefClose(0)));
    EXPECT_EQ(2,(matrix.getNumRefClose(5)));
    EXPECT_EQ(2,(matrix.getNumRefClose(10)));
    EXPECT_EQ(3,(matrix.getNumRefClose(15)));
    EXPECT_EQ(2,(matrix.getNumRefClose(20)));
}

TEST(TestOptiRefMatrix, isCloseFit) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
    
    bool isFit;
    vector<long long> fitSeqs = matrix.getFitSeqs();
    string Expected_ReturnResults = "024689101213141721232426293036394243464950525457";
    string ReturnResults = "";
    for (long long i = 0; i < fitSeqs.size(); i++) { ReturnResults += toString(fitSeqs[i]);  }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    //check closeness
    EXPECT_EQ(true,(matrix.isClose(0, 8)));
    EXPECT_EQ(true,(matrix.isClose(1, 28)));
    EXPECT_EQ(true,(matrix.isClose(2, 44)));
    EXPECT_EQ(true,(matrix.isClose(15, 42)));
    EXPECT_EQ(true,(matrix.isClose(35, 36)));
    
    //check not close
    EXPECT_EQ(false,(matrix.isClose(57, 8)));
    EXPECT_EQ(false,(matrix.isClose(47, 28)));
    EXPECT_EQ(false,(matrix.isClose(32, 44)));
    EXPECT_EQ(false,(matrix.isClose(23, 42)));
    EXPECT_EQ(false,(matrix.isClose(12, 36)));
    
    //assumes first value is a fitSeq
    EXPECT_EQ(true,(matrix.isCloseFit(0, 8, isFit))); //both fit and close
    EXPECT_EQ(true, isFit);
    EXPECT_EQ(false,(matrix.isCloseFit(2, 28, isFit))); //not fit
    EXPECT_EQ(false, isFit);
    EXPECT_EQ(false,(matrix.isCloseFit(3, 20, isFit))); //not fit, but close
    EXPECT_EQ(false, isFit);
    EXPECT_EQ(false,(matrix.isCloseFit(13, 42, isFit))); //both fit, not close
    EXPECT_EQ(true, isFit);
    EXPECT_EQ(false,(matrix.isCloseFit(30, 36, isFit))); //both fit, not close
    EXPECT_EQ(true, isFit);
}


TEST(TestOptiRefMatrix, getCloseFitSeqs) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
 
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
 
    //"024689101213141721232426293036394243464950525457";
 
    //17	GQY1XT001BJ4H6,..,GQY1XT001CW8RQ	11	32	52	55	57
    string Expected_ReturnResults = ""; Expected_ReturnResults += "52"; Expected_ReturnResults += "57";
    set<long long> temp = matrix.getCloseFitSeqs(17);
    string ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    //50	GQY1XT001EK1FO	13
    Expected_ReturnResults = ""; Expected_ReturnResults += "13";
    temp = matrix.getCloseFitSeqs(50);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    //52	GQY1XT001ENMKV	11	17	32	55	57
    Expected_ReturnResults = ""; Expected_ReturnResults += "17"; Expected_ReturnResults += "57";
    temp = matrix.getCloseFitSeqs(52);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    //36	GQY1XT001DHDV0,GQY1XT001B0UFF	15	24	35	38
    Expected_ReturnResults = ""; Expected_ReturnResults += "24";
    temp = matrix.getCloseFitSeqs(36);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    //21	GQY1XT001BUMO0	26	42	46
    Expected_ReturnResults = ""; Expected_ReturnResults += "26"; Expected_ReturnResults += "42"; Expected_ReturnResults += "46";
    temp = matrix.getCloseFitSeqs(21);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
}

TEST(TestOptiRefMatrix, extractRefMatrix) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
 
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
    OptiData* refMatrix = matrix.extractRefMatrix();
    
    EXPECT_EQ(54,(refMatrix->print(cout)));
    EXPECT_EQ(54,(refMatrix->getNumDists()));
}

TEST(TestOptiRefMatrix, extractMatrixSubset) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
    vector<long long> temp = matrix.getFitSeqs();
    set<long long> fitSeqs;
    for (long long i = 0; i < temp.size(); i++) { fitSeqs.insert(temp[i]); }
    OptiData* fitMatrix = matrix.extractMatrixSubset(fitSeqs);
    
    EXPECT_EQ(34,(fitMatrix->print(cout)));
    EXPECT_EQ(34,(fitMatrix->getNumDists()));
}


TEST(TestOptiRefMatrix, getFitListSingle) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);
    
    //maps names to index in closeness matrix
    ListVector* fitListSingle = matrix.getFitListSingle();
    
    //check bin 0
    string bin = fitListSingle->get(0);
    EXPECT_EQ("GQY1XT001AD34Z", bin);
    
    bin = fitListSingle->get(10);
    EXPECT_EQ("GQY1XT001CH9UX,GQY1XT001C80OT,GQY1XT001BEIF2,GQY1XT001DFU9M,GQY1XT001DNJRS", bin);
    
    bin = fitListSingle->get(16);
    EXPECT_EQ("GQY1XT001EACH9", bin);
    
    bin = fitListSingle->get(3);
    EXPECT_EQ("GQY1XT001B8C4W,GQY1XT001DBTGA,GQY1XT001B4VQ6", bin);
    
    bin = fitListSingle->get(7);
    EXPECT_EQ("GQY1XT001C4UVG", bin);
}

TEST(TestOptiRefMatrix, randomizeRefs) {
    MothurOut* m; m = MothurOut::getInstance();
    m->setRandomSeed(123456); //stabilize radomization
    
    TestOptiRefMatrix testOMatrix;
    OptiRefMatrix matrix(testOMatrix.columnFile, "column", testOMatrix.filenames[1], "name", 0.03, 50);

    matrix.randomizeRefs();
    
    vector<long long> refSeqs = matrix.getRefSeqs();
    string Expected_ReturnResults = "0135678112023242629313236394041434547525356";
    string ReturnResults = "";
    for (long long i = 0; i < refSeqs.size(); i++) { ReturnResults += toString(refSeqs[i]);  }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    matrix.randomizeRefs();
    
    refSeqs = matrix.getRefSeqs();
    Expected_ReturnResults = "3467891315161718202122232426273132333536373841474950";
    ReturnResults = "";
    for (long long i = 0; i < refSeqs.size(); i++) { ReturnResults += toString(refSeqs[i]);  }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
}

/**************************************************************************************************/


