#ifndef MAIN_TEST
#define MAIN_TEST

//
//  main.cpp
//  TestMothur
//
//  Created by Sarah Westcott on 3/23/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//


#include "mothurout.h"
#include "currentfile.h"
#include "commandfactory.hpp"
#include "gtest/gtest.h"

/* Test Naming Structure:
 
 Test_OptionalSubCategory_TestClass
 
 Test_Container_Sequence
 Test_Calcs_ClusterCalcs
 
 Makes it easy to filter tests
 
 ::testing::GTEST_FLAG(filter) = "Test_Container_*";
 ::testing::GTEST_FLAG(filter) = "Test_*";
 ::testing::GTEST_FLAG(filter) = "Test_Calcs*";
 
 Test_TrimOligos
 
 
*/

CommandFactory* CommandFactory::_uniqueInstance;
CurrentFile* CurrentFile::instance;
MothurOut* MothurOut::_uniqueInstance;

int main(int argc, char **argv) {
    MothurOut* m; m = MothurOut::getInstance();
    CurrentFile* current; current = CurrentFile::getInstance();
    
    current->setTestFilePath("/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/");
    string pathname = current->getProgramPath();
    if (pathname != "") {
        //add / to name if needed
        string lastChar = pathname.substr(pathname.length()-1);
        if (lastChar != PATH_SEPARATOR) { pathname += PATH_SEPARATOR; }
    }
    
    current->setTestFilePath(pathname);
    
    ::testing::GTEST_FLAG(filter) = "Test_TrimOligos*";
    ::testing::InitGoogleTest(&argc, argv);
    
    int value = RUN_ALL_TESTS();
    return value;
}

#endif
