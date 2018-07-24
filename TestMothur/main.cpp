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
#include "gtest/gtest.h"


#define UNIT_TEST

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
    
    ::testing::InitGoogleTest(&argc, argv);
    
    return RUN_ALL_TESTS();
}

#endif
