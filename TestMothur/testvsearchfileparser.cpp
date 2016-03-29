//
//  testvsearchfileparser.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/24/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "catch.hpp"
#include "testvsearchfileparser.h"

/**************************************************************************************************/
TestVsearchFileParser::TestVsearchFileParser() {  //setup
    m = MothurOut::getInstance();
    TestDataSet data;
    filenames = data.getSubsetFNGFiles(10);
    ct = data.getCountTable();
}
/**************************************************************************************************/
TestVsearchFileParser::~TestVsearchFileParser() {
    delete ct;
    for (int i = 0; i < filenames.size(); i++) { m->mothurRemove(filenames[i]); } //teardown
}
/**************************************************************************************************/

TEST_CASE("Testing VsearchParser Class") {
    TestVsearchFileParser testVParser;
    VsearchFileParser vsearchParser(testVParser.filenames[0], testVParser.filenames[1], "name");
    
    SECTION("CreateVsearchFasta") {
        INFO("Using First 10 sequences of final.fasta and final.names") // Only appears on a FAIL
        
        CAPTURE(vsearchParser.getVsearchFile()); // Displays this variable on a FAIL
        
        CHECK(vsearchParser.getVsearchFile() == "tempSeqs.sorted.fasta.temp");
        
        ifstream in;
        testVParser.m->openInputFile(vsearchParser.getVsearchFile(), in);
        
        while (!in.eof()) {
            Sequence seq(in); testVParser.m->gobble(in);
            
            vector<string> pieces;
            string name = seq.getName();
            testVParser.m->splitAtChar(name, pieces, '=');
            string abundString = pieces[1].substr(0, pieces[1].length()-1);
            int abund = 0;
            testVParser.m->mothurConvert(abundString, abund);
            int totalSeqs = testVParser.ct->getNumSeqs(seq.getName());
            
            CHECK(abund == totalSeqs);
        }
    }
    
    SECTION("Remove Abundances") {
        INFO("Using GQY1XT001C44N8/size=3677/") // Only appears on a FAIL
        string seqName = "GQY1XT001C44N8/size=3677/";
        CAPTURE(testVParser.removeAbundances(seqName)); // Displays this variable on a FAIL
        
        CHECK(testVParser.removeAbundances(seqName) == "GQY1XT001C44N8");
    }
}
/**************************************************************************************************/
