//
//  testbiominfocommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/18/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//

#include "testbiominfocommand.h"
#include "catch.hpp"


TEST_CASE("Testing Biom.Info Command") {
    TestBiomInfoCommand tBiom;
    string input = "";
    
    SECTION("Test getName function") {
        INFO("Using id:B") // Only appears on a FAIL
        string input = "id:B";
        CAPTURE(tBiom.getName(input)); // Displays this variable on a FAIL
        
        CHECK(tBiom.getName(input) == "B");
    }
    
    SECTION("Test getTaxonomy function") {
        INFO("Using taxonomy:k__Bacteria,p__Firmicutes,c__Bacilli,o__Turicibacterales,f__Turicibacteraceae,g__Turicibacter,s__	bootstrap:100,100,100,100,100,100,100") // Only appears on a FAIL
        string tax = "taxonomy:k__Bacteria,p__Firmicutes,c__Bacilli,o__Turicibacterales,f__Turicibacteraceae,g__Turicibacter,s__";
        string boot = "bootstrap:100,100,100,100,100,100,100";
        CAPTURE(tBiom.getTaxonomy(tax, boot)); // Displays this variable on a FAIL
        
        CHECK(tBiom.getTaxonomy(tax, boot) == "k__Bacteria(100);p__Firmicutes(100);c__Bacilli(100);o__Turicibacterales(100);f__Turicibacteraceae(100);g__Turicibacter(100);s__(100);");
    }
    
    SECTION("Test readRows") {
        INFO("Using columns:[{id:A, metadata:null},{id:B, metadata:null},{id:C, metadata:null}]") // Only appears on a FAIL
        string input = "columns:[{id:A, metadata:null},{id:B, metadata:null},{id:C, metadata:null}]";
        int numRows = 0;
        bool hasTaxonomy = true;
        //CAPTURE(tBiom.readRows(input, numRows, hasTaxonomy)); // Displays this variable on a FAIL
        
        CHECK(tBiom.readRows(input, numRows, hasTaxonomy)[0][0] == "A");
        CHECK(numRows == 3);
        CHECK(hasTaxonomy == false);
    }
    
    SECTION("Test getDims") {
        INFO("Using line, shape: [28,3]")
        string input = "shape: [28,3]";
        int numRows = 0;
        int numCols = 0;
        CAPTURE(tBiom.getDims(input, numRows, numCols)); // Displays this variable on a FAIL
        
        CHECK(tBiom.getDims(input, numRows, numCols) == 0);
        CHECK(numRows == 28);
        CHECK(numCols == 3);
    }
        
    SECTION("Test readData") {
        INFO("Using line, data:  [[0,2,5], [1,2,5], [2,1,2], [3,1,1], [4,1,1], [5,0,18], [5,1,12], [6,0,15], [6,1,4], [7,0,1], [7,1,1], [8,1,1], [9,0,2], [9,1,6], [9,2,4], [10,1,2], [11,0,5], [11,1,1], [11,2,4], [12,0,1], [13,0,1], [13,1,2], [14,1,2], [15,1,5], [16,0,8], [16,1,1], [17,1,2], [18,0,13], [19,0,2], [19,1,1], [20,0,15], [20,1,27], [20,2,11], [21,1,10], [22,2,18], [23,2,5], [24,0,1], [24,2,20], [25,1,2], [25,2,2], [26,0,1], [26,1,1], [27,0,1]]")
        string input = "data:  [[0,2,5], [1,2,5], [2,1,2], [3,1,1], [4,1,1], [5,0,18], [5,1,12], [6,0,15], [6,1,4], [7,0,1], [7,1,1], [8,1,1], [9,0,2], [9,1,6], [9,2,4], [10,1,2], [11,0,5], [11,1,1], [11,2,4], [12,0,1], [13,0,1], [13,1,2], [14,1,2], [15,1,5], [16,0,8], [16,1,1], [17,1,2], [18,0,13], [19,0,2], [19,1,1], [20,0,15], [20,1,27], [20,2,11], [21,1,10], [22,2,18], [23,2,5], [24,0,1], [24,2,20], [25,1,2], [25,2,2], [26,0,1], [26,1,1], [27,0,1]]""";
        string matrixFormat = "sparse";
        string matrixElementType = "int";
        vector<string> groupNames; groupNames.push_back("A"); groupNames.push_back("B"); groupNames.push_back("C");
        int numOtus = 28;
        //CAPTURE(tBiom.readData(matrixFormat, input, matrixElementType, groupNames, numOtus)); // Displays this variable on a FAIL
        int abund = tBiom.readData(matrixFormat, input, matrixElementType, groupNames, numOtus)[0]->getAbundance(9);
        CHECK(abund == 2);
    }
    
    
    //more tests need to be added - just a start to set up testing project and model
}

