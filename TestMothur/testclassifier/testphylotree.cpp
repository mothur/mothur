//
//  testphylotree.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/29/16.
//  Copyright © 2016 Schloss Lab. All rights reserved.
//

#include "catch.hpp"
#include "testphylotree.hpp"

/**************************************************************************************************/
TestPhyloTree::TestPhyloTree() {  //setup
    m = MothurOut::getInstance();
    
    string tax1WithSpaces = "Bacteria(100);Bacteroidetes 7(100);Bacteroidia(100);Bacteroidales(100);S24-7(100);";
    string tax2WithSpaces = "Bacteria(100);Bacteroidetes 7(100);Bacteroidia(98);Bacteroidales(98);Bacteroidaceae(98);Bacteroides(98);";
    string tax3WithSpaces = "Bacteria(100);Firmicutes(100);Clostridia B(100);Clostridiales(100);Lachnospiraceae(100);Blautia(92);";
    string tax4WithSpaces = "Bacteria(100);Firmicutes(100);Clostridia B(100);Clostridiales(100);Ruminococcaceae(100);Anaerotruncus(100);";
    string tax5WithSpaces = "Bacteria(100);Firmicutes(100);Clostridia B(100);Clostridiales(100);Lachnospiraceae(100);Incertae_Sedis(97);";
    
    string tax1WithOutSpaces = "Bacteria(100);Bacteroidetes(100);Bacteroidia(100);Bacteroidales(100);S24-7(100);";
    string tax2WithOutSpaces = "Bacteria(100);Bacteroidetes(100);Bacteroidia(98);Bacteroidales(98);Bacteroidaceae(98);Bacteroides(98);";
    string tax3WithOutSpaces = "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Blautia(92);";
    string tax4WithOutSpaces = "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(100);Anaerotruncus(100);";
    string tax5WithOutSpaces = "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Incertae_Sedis(97);";
    
    phylo.addSeqToTree("seq1", tax1WithSpaces);
    phylo.addSeqToTree("seq2", tax2WithSpaces);
    phylo.addSeqToTree("seq3", tax3WithSpaces);
    phylo.addSeqToTree("seq4", tax4WithSpaces);
    phylo.addSeqToTree("seq5", tax5WithSpaces);
    phylo.addSeqToTree("seq6", tax1WithOutSpaces);
    phylo.addSeqToTree("seq7", tax2WithOutSpaces);
    phylo.addSeqToTree("seq8", tax3WithOutSpaces);
    phylo.addSeqToTree("seq9", tax4WithOutSpaces);
    phylo.addSeqToTree("seq10", tax5WithOutSpaces);
}
/**************************************************************************************************/
TestPhyloTree::~TestPhyloTree() {}
/**************************************************************************************************/

TEST_CASE("Testing PhyloTree Class") {
    TestPhyloTree testPTree;
    
    string tax1WithSpaces = "Bacteria(100);Bacteroidetes 7(100);Bacteroidia(100);Bacteroidales(100);S24-7(100);";
    string tax2WithSpaces = "Bacteria(100);Bacteroidetes 7(100);Bacteroidia(98);Bacteroidales(98);Bacteroidaceae(98);Bacteroides(98);";
    string tax3WithSpaces = "Bacteria(100);Firmicutes(100);Clostridia B(100);Clostridiales(100);Lachnospiraceae(100);Blautia(92);";
    string tax4WithSpaces = "Bacteria(100);Firmicutes(100);Clostridia B(100);Clostridiales(100);Ruminococcaceae(100);Anaerotruncus(100);";
    string tax5WithSpaces = "Bacteria(100);Firmicutes(100);Clostridia B(100);Clostridiales(100);Lachnospiraceae(100);Incertae_Sedis(97);";
    
    string tax1WithOutSpaces = "Bacteria(100);Bacteroidetes(100);Bacteroidia(100);Bacteroidales(100);S24-7(100);";
    string tax2WithOutSpaces = "Bacteria(100);Bacteroidetes(100);Bacteroidia(98);Bacteroidales(98);Bacteroidaceae(98);Bacteroides(98);";
    string tax3WithOutSpaces = "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Blautia(92);";
    string tax4WithOutSpaces = "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(100);Anaerotruncus(100);";
    string tax5WithOutSpaces = "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Incertae_Sedis(97);";
    
    SECTION("Add Sequences to Tree") {
        INFO("Using taxonomies with and without spaces") // Only appears on a FAIL
        
        CAPTURE(testPTree.addSeqToTree("seq1", tax1WithSpaces));
        CHECK(testPTree.addSeqToTree("seq1", tax1WithSpaces) == 5);
        
        CAPTURE(testPTree.addSeqToTree("seq2", tax2WithSpaces));
        CHECK(testPTree.addSeqToTree("seq2", tax2WithSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq3", tax3WithSpaces));
        CHECK(testPTree.addSeqToTree("seq3", tax3WithSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq4", tax4WithSpaces));
        CHECK(testPTree.addSeqToTree("seq4", tax4WithSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq5", tax5WithSpaces));
        CHECK(testPTree.addSeqToTree("seq5", tax5WithSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq6", tax1WithOutSpaces));
        CHECK(testPTree.addSeqToTree("seq6", tax1WithOutSpaces) == 5);
        
        CAPTURE(testPTree.addSeqToTree("seq7", tax2WithOutSpaces));
        CHECK(testPTree.addSeqToTree("seq7", tax2WithOutSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq8", tax3WithOutSpaces));
        CHECK(testPTree.addSeqToTree("seq8", tax3WithOutSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq9", tax4WithOutSpaces));
        CHECK(testPTree.addSeqToTree("seq9", tax4WithOutSpaces) == 6);
        
        CAPTURE(testPTree.addSeqToTree("seq10", tax5WithOutSpaces));
        CHECK(testPTree.addSeqToTree("seq10", tax5WithOutSpaces) == 6);
    }
    
    SECTION("Get Seqs") {
        INFO("Using taxonomies with and without spaces") // Only appears on a FAIL
               
        CAPTURE(testPTree.phylo.getSeqs("Bacteroidetes 7").size());
        CHECK((testPTree.phylo.getSeqs("Bacteroidetes 7").size()) == 2);
        
        vector<string> Bacteroidetes_7 = testPTree.phylo.getSeqs("Bacteroidetes 7");
        CHECK(Bacteroidetes_7[0] == "seq1");
        CHECK(Bacteroidetes_7[1] == "seq2");
        
        CAPTURE(testPTree.phylo.getSeqs("Clostridia").size());
        CHECK((testPTree.phylo.getSeqs("Clostridia").size()) == 3);
        
        vector<string> Clostridia = testPTree.phylo.getSeqs("Clostridia");
        CHECK(Clostridia[0] == "seq8");
        CHECK(Clostridia[1] == "seq9");
        CHECK(Clostridia[2] == "seq10");
    }
   
    SECTION("Get Genus Totals") {
        INFO("Using taxonomies with and without spaces") // Only appears on a FAIL
        
        CAPTURE(testPTree.phylo.getGenusTotals().size());
        CHECK(testPTree.phylo.getGenusTotals().size() == 10);
    }

    SECTION("Get Full Taxonomy") {
        INFO("Using taxonomies with and without spaces") // Only appears on a FAIL
        
        CAPTURE(testPTree.phylo.getFullTaxonomy("seq1"));
        CHECK(testPTree.phylo.getFullTaxonomy("seq1") == "Bacteria;Bacteroidetes 7;Bacteroidia;Bacteroidales;S24-7;");
        
        CAPTURE(testPTree.phylo.getFullTaxonomy("seq10"));
        CHECK(testPTree.phylo.getFullTaxonomy("seq1") == "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Incertae_Sedis;");
    }
    
}
/**************************************************************************************************/
