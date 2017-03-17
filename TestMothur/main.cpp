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
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
