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

TEST_CASE("Testing VsearchParser Class") {
    MothurOut* m = MothurOut::getInstance();
    TestDataSet data;
    vector<string> filenames = data.getSubsetFNGFiles(10);
    CountTable* ct = data.getCountTable();
    
    VsearchFileParser vsearchParser(filenames[0], filenames[1], "name");
    TestVsearchFileParser testVParser;
    
    SECTION("CreateVsearchFasta") {
        INFO("Using First 10 sequences of final.fasta and final.names") // Only appears on a FAIL
        
        CAPTURE(vsearchParser.getVsearchFile()); // Displays this variable on a FAIL
        
        CHECK(vsearchParser.getVsearchFile() == "tempSeqs.sorted.fasta.temp");
        
        ifstream in;
        m->openInputFile(vsearchParser.getVsearchFile(), in);
        
        while (!in.eof()) {
            Sequence seq(in); m->gobble(in);
            
            vector<string> pieces;
            string name = seq.getName();
            m->splitAtChar(name, pieces, '=');
            string abundString = pieces[1].substr(0, pieces[1].length()-1);
            int abund = 0;
            m->mothurConvert(abundString, abund);
            int totalSeqs = ct->getNumSeqs(seq.getName());
            
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
