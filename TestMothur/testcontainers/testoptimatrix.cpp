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
    CurrentFile* current; current = CurrentFile::getInstance();
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(); //Fasta, name, group returned
    
    columnFile = data.getSubsetFNGDistFile();
    phylipFile = data.getSubsetFNGPhylipDistFile();
    
}
/**************************************************************************************************/
TestOptiMatrix::~TestOptiMatrix() {}
/**************************************************************************************************/
////distfile, dupsFile, dupsFormat, distFormat, cutoff, sim
TEST(TestOptiMatrix, readColumn) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", 0.03, false);
    
    EXPECT_EQ(160,(matrix.print(cout)));
    EXPECT_EQ(160,(matrix.getNumDists()));
}

TEST(TestOptiMatrix, readPhylip) {
    TestOptiMatrix testOMatrix;
    OptiMatrix pmatrix(testOMatrix.phylipFile, "", "", "phylip", 0.03, false);

    EXPECT_EQ(160,(pmatrix.print(cout)));
    EXPECT_EQ(160,(pmatrix.getNumDists()));
    
}

TEST(TestOptiMatrix, getNumCLose) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", 0.03, false);
    
    EXPECT_EQ(1,(matrix.getNumClose(0)));
    EXPECT_EQ(2,(matrix.getNumClose(5)));
    EXPECT_EQ(3,(matrix.getNumClose(10)));
    EXPECT_EQ(7,(matrix.getNumClose(15)));
    EXPECT_EQ(2,(matrix.getNumClose(20)));
}

TEST(TestOptiMatrix, isClose) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", 0.03, false);
    
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
}

TEST(TestOptiMatrix, getCloseSeqs) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", 0.03, false);
    
    //11	GQY1XT001B04KZ,GQY1XT001EBRFH	17	32	52	55	57
    string Expected_ReturnResults = ""; Expected_ReturnResults += "17"; Expected_ReturnResults += "32"; Expected_ReturnResults += "52"; Expected_ReturnResults += "55"; Expected_ReturnResults += "57";
    set<long long> temp = matrix.getCloseSeqs(11);
    string ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
    
    //21	GQY1XT001BUMO0	26	42	46
    Expected_ReturnResults = ""; Expected_ReturnResults += "26"; Expected_ReturnResults += "42"; Expected_ReturnResults += "46";
    temp = matrix.getCloseSeqs(21);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);

    //31	GQY1XT001CVCKG,GQY1XT001BO8Z9	20	27
    Expected_ReturnResults = ""; Expected_ReturnResults += "20"; Expected_ReturnResults += "27";
    temp = matrix.getCloseSeqs(31);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);

    //41	GQY1XT001DY3E7	19	29
    Expected_ReturnResults = ""; Expected_ReturnResults += "19"; Expected_ReturnResults += "29";
    temp = matrix.getCloseSeqs(41);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);

    //51	GQY1XT001EN363,GQY1XT001B0ZKY,GQY1XT001BCPXE,GQY1XT001BEKE1,GQY1XT001D25E1,GQY1XT001EWORZ,GQY1XT001AQB9P,GQY1XT001CEFI4	49
    Expected_ReturnResults = ""; Expected_ReturnResults += "49";
    temp = matrix.getCloseSeqs(51);
    ReturnResults = "";
    for (set<long long>::iterator it = temp.begin(); it != temp.end(); it++) { ReturnResults += toString(*it); }
    
    EXPECT_EQ(Expected_ReturnResults, ReturnResults);
}

TEST(TestOptiMatrix, getNameIndexMap) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", 0.03, false);
    
    //maps names to index in closeness matrix
    map<string, long long> nameIndexMap = matrix.getNameIndexMap();
    
    //check nameMap
    EXPECT_EQ(0,nameIndexMap["GQY1XT001A4DGI"]);
    EXPECT_EQ(39,nameIndexMap["GQY1XT001DRYVA"]);
    EXPECT_EQ(44,nameIndexMap["GQY1XT001E23UK"]);
    EXPECT_EQ(52,nameIndexMap["GQY1XT001ENMKV"]);
    EXPECT_EQ(48,nameIndexMap["GQY1XT001EJAUJ"]);
    
    EXPECT_EQ("GQY1XT001A4DGI",(matrix.getName(0)));
    EXPECT_EQ("GQY1XT001DRYVA",(matrix.getName(39)));
    EXPECT_EQ("GQY1XT001E23UK",(matrix.getName(44)));
    EXPECT_EQ("GQY1XT001ENMKV",(matrix.getName(52)));
    EXPECT_EQ("GQY1XT001EJAUJ",(matrix.getName(48)));
}

TEST(TestOptiMatrix, getListSingle) {
    TestOptiMatrix testOMatrix;
    OptiMatrix matrix(testOMatrix.columnFile, testOMatrix.filenames[1], "name", "column", 0.03, false);
    
    //maps names to index in closeness matrix
    ListVector* listSingle = matrix.getListSingle();
    
    //check bin 0
    string bin = listSingle->get(0);
    EXPECT_EQ("GQY1XT001AD34Z", bin);
    
    bin = listSingle->get(10);
    EXPECT_EQ("GQY1XT001BRLCO", bin);

    bin = listSingle->get(18);
    EXPECT_EQ("GQY1XT001CKAUI", bin);

    bin = listSingle->get(3);
    EXPECT_EQ("GQY1XT001AOSH9,GQY1XT001BLJ4I,GQY1XT001BNIJQ,GQY1XT001CT9JB,GQY1XT001DCPGQ,GQY1XT001DY88Y,GQY1XT001AHO0L,GQY1XT001DRMZK,GQY1XT001DIXY7,GQY1XT001CDBZ1,GQY1XT001B8C47,GQY1XT001A71WZ,GQY1XT001D41QJ,GQY1XT001BAMTS", bin);

    bin = listSingle->get(7);
    EXPECT_EQ("GQY1XT001B8UKY", bin);
}
/**************************************************************************************************/


