//
//  testphylotree.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/29/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//

#include "testphylotree.hpp"

/**************************************************************************************************/
TestPhyloTree::TestPhyloTree() {  //setup
    m = MothurOut::getInstance();
}
/**************************************************************************************************/
TestPhyloTree::~TestPhyloTree() {}
/**************************************************************************************************/

TEST_CASE("Testing PhyloTree Class") {
    TestPhyloTree testPTree;
    PhyloTree phylo;
    
    SECTION("Add Sequences to Tree") {
        INFO("Using taxonomies with and without spaces") // Only appears on a FAIL
        
        CAPTURE(vsearchParser.getVsearchFile()); // Displays this variable on a FAIL
        
        CHECK(vsearchParser.getVsearchFile() == "tempSeqs.txt.sorted.fasta.temp");
        
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
            int totalSeqs = testVParser.ct->getNumSeqs(testVParser.removeAbundances(name));
            
            CHECK(abund == totalSeqs);
        }
        in.close();
        testVParser.m->mothurRemove("tempSeqs.txt.sorted.fasta.temp");
    }
    
}
/**************************************************************************************************/
