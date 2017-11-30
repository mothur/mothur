//
//  testrenamefilecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/4/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "testrenamefilecommand.h"
#include "dataset.h"


/**************************************************************************************************/
TestRenameFileCommand::TestRenameFileCommand() {  //setup
    m = MothurOut::getInstance();
    TestDataSet data;
    filenames = data.getSubsetFNGFiles();
}
/**************************************************************************************************/
TestRenameFileCommand::~TestRenameFileCommand() {}
/**************************************************************************************************
TEST_CASE("Testing RenameFileCommand Class") {
    TestRenameFileCommand testRename;
    
    SECTION("Testing GetNewName - with prefix") {
        INFO("Using prefix=greatData") // Only appears on a FAIL
        
        testRename.prefix = "greatData";
        testRename.mothurGenerated = true;
        
        CAPTURE(testRename.getNewName(testRename.filenames[0], "fasta")); // Displays this variable on a FAIL
        
        CHECK(testRename.getNewName(testRename.filenames[0], "fasta") == "greatData.txt");
        
        testRename.filenames[0] = testRename.getNewName(testRename.filenames[0], "fasta"); //for teardown
    }
    
    SECTION("Testing GetNewName - with user name") {
        INFO("Using prefix=greatData") // Only appears on a FAIL
        
        testRename.outputfile = "greatData.fasta";
        testRename.mothurGenerated = false;
        
        CAPTURE(testRename.getNewName(testRename.filenames[0], "fasta")); // Displays this variable on a FAIL
        
        CHECK(testRename.getNewName(testRename.filenames[0], "fasta") == "greatData.fasta");
        
        testRename.filenames[0] = testRename.getNewName(testRename.filenames[0], "fasta");  //for teardown
    }
    
    
    SECTION("Testing RenameOrCopy - deleteOld=false") {
        INFO("Uses mothur rename function to move or system command to copy.") // Only appears on a FAIL
        
        testRename.deleteOld = false;
        
        testRename.renameOrCopy(testRename.filenames[0], "greatData.new.fasta");
        
        ifstream in, in2;
        bool ableToOpen = testRename.util.openInputFile("greatData.new.fasta", in);
        in.close();
        
        CAPTURE(ableToOpen);
            
        CHECK(ableToOpen == 0);
        
        bool ableToOpen2 = testRename.util.openInputFile(testRename.filenames[0], in2);
        in2.close();
        
        CAPTURE(ableToOpen2);
        
        CHECK(ableToOpen2 == 0);
        
        testRename.util.mothurRemove("greatData.new.fasta");
    }
    
    SECTION("Testing RenameOrCopy - deleteOld=true") {
        INFO("Uses mothur rename function to move or system command to copy.") // Only appears on a FAIL
        
        testRename.deleteOld = true;
        
        testRename.renameOrCopy(testRename.filenames[0], "greatData.new.fasta");
        
        ifstream in, in2;
        bool ableToOpen = testRename.util.openInputFile("greatData.new.fasta", in);
        in.close();
        
        CAPTURE(ableToOpen);
        
        CHECK(ableToOpen == 0);
        
        testRename.filenames[0] = testRename.getNewName(testRename.filenames[0], "fasta");  //for teardown
    }
}*/
/**************************************************************************************************/
