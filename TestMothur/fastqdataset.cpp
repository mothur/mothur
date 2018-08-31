//
//  fastqdataset.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/31/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "fastqdataset.h"

/***********************************************************************/
TestFastqDataSet::TestFastqDataSet() {  m = MothurOut::getInstance();  }
/***********************************************************************/
vector<string> TestFastqDataSet::getSubsetFRFastq(int numSeqs) {
    
    fillForwardFastq();
    fillReverseFastq();
    ofstream out, out2;
    util.openOutputFile("tempForward.txt", out); util.openOutputFile("tempReverse.txt", out2);
    for (int i = 0; i < numSeqs; i++) {
        ffastqReads[i].printFastq(out);
        rfastqReads[i].printFastq(out2);
    }
    
    vector<string> filenames; filenames.push_back("tempForward.txt"); filenames.push_back("tempReverse.txt");
    
    return filenames;
}
/***********************************************************************/
void TestFastqDataSet::fillForwardFastq() {
    ffastqReads.clear();
    
    //read info from stable file
    string testfile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/F8D0_S345_L001_R1_001.fastq";
    
    ifstream in;
    util.openInputFile(testfile, in);
    
    int count = 0; bool ignore = false; string format = "illumina1.8+";
    while (!in.eof()) {
        if (m->getControl_pressed()) { break; }
        
        if (count < 2000) {
            FastqRead read(in, ignore, format); util.gobble(in);
            if (!ignore) { ffastqReads.push_back(read);  count++; }
        }else { break; }
        
    }
    in.close();
}
/***********************************************************************/
void TestFastqDataSet::fillReverseFastq() {
    rfastqReads.clear();
    
    //read info from stable file
    string testfile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/F8D0_S345_L001_R2_001.fastq";
    
    ifstream in;
    util.openInputFile(testfile, in);
    
    int count = 0; bool ignore = false; string format = "illumina1.8+";
    while (!in.eof()) {
        if (m->getControl_pressed()) { break; }
        
        if (count < 2000) {
            FastqRead read(in, ignore, format); util.gobble(in);
            if (!ignore) { ffastqReads.push_back(read);  count++; }
        }else { break; }
        
    }
    in.close();
}
/***********************************************************************/
