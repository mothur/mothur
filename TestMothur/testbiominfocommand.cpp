//
//  testbiominfocommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/18/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "testbiominfocommand.h"

/**************************************************************************************************/
TestBiomInfoCommand::TestBiomInfoCommand() {
    m = MothurOut::getInstance();
}
/**************************************************************************************************/
TestBiomInfoCommand::~TestBiomInfoCommand() { }
/**************************************************************************************************/

TEST(Test_Command_BiomInfo, getDims) {
    TestBiomInfoCommand testTrim;
    
    string input = "shape: [28,3]";
    int numRows = 0; int numCols = 0;
    testTrim.getDims(input, numRows, numCols);
    
    EXPECT_EQ(28, numRows);
    EXPECT_EQ(3, numCols);
}

TEST(Test_Command_BiomInfo, getName) {
    TestBiomInfoCommand testTrim;
    
    string input = "id:B";
    EXPECT_EQ("B", testTrim.getName(input));
}

TEST(Test_Command_BiomInfo, getTaxonomy) {
    TestBiomInfoCommand testTrim;
    
    string tax = "taxonomy:k__Bacteria,p__Firmicutes,c__Bacilli,o__Turicibacterales,f__Turicibacteraceae,g__Turicibacter,s__";
    string boot = "bootstrap:100,100,100,100,100,100,100";
    EXPECT_EQ("k__Bacteria(100);p__Firmicutes(100);c__Bacilli(100);o__Turicibacterales(100);f__Turicibacteraceae(100);g__Turicibacter(100);s__(100);", testTrim.getTaxonomy(tax, boot));
}

TEST(Test_Command_BiomInfo, readRows) {
    TestBiomInfoCommand testTrim;
    
    string input = "columns:[{id:A, metadata:null},{id:B, metadata:null},{id:C, metadata:null}]";
    int numRows = 0;
    bool hasTaxonomy = true;
    EXPECT_EQ("A", testTrim.readRows(input, numRows, hasTaxonomy)[0][0]);
    EXPECT_EQ(3, numRows);
    EXPECT_EQ(false, hasTaxonomy);
}

TEST(Test_Command_BiomInfo, readData) {
    TestBiomInfoCommand testTrim;
    
    string input = "data:  [[0,2,5], [1,2,5], [2,1,2], [3,1,1], [4,1,1], [5,0,18], [5,1,12], [6,0,15], [6,1,4], [7,0,1], [7,1,1], [8,1,1], [9,0,2], [9,1,6], [9,2,4], [10,1,2], [11,0,5], [11,1,1], [11,2,4], [12,0,1], [13,0,1], [13,1,2], [14,1,2], [15,1,5], [16,0,8], [16,1,1], [17,1,2], [18,0,13], [19,0,2], [19,1,1], [20,0,15], [20,1,27], [20,2,11], [21,1,10], [22,2,18], [23,2,5], [24,0,1], [24,2,20], [25,1,2], [25,2,2], [26,0,1], [26,1,1], [27,0,1]]""";
    string matrixFormat = "sparse";
    string matrixElementType = "int";
    vector<string> groupNames; groupNames.push_back("A"); groupNames.push_back("B"); groupNames.push_back("C");
    int numOtus = 28;
    
    EXPECT_EQ(12, testTrim.readData(matrixFormat, input, matrixElementType, groupNames, numOtus)->getOTUTotal(9));
}

/**************************************************************************************************/


