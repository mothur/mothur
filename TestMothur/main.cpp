//
//  main.cpp
//  TestMothur
//
//  Created by Sarah Westcott on 3/23/15.
//  Copyright (c) 2015 Schloss Lab. All rights reserved.
//


#include "mothurout.h"
#include "gtest/gtest.h"

#define UNIT_TEST

int main(int argc, char **argv) {
    MothurOut* m; m = MothurOut::getInstance();
    string pathname = m->mothurProgramPath;
    if (pathname != "") {
        //add / to name if needed
        string lastChar = pathname.substr(pathname.length()-1);
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        if (lastChar != "/") { pathname += "/TestMothur/TestFiles/"; }
#else
        if (lastChar != "\\") { pathname += "\\TestMothur\\TestFiles\\"; }
#endif
    }
    
    m->setTestFilePath(pathname);
    
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
