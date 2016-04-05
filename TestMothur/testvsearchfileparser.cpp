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
    filenames = data.getSubsetFNGFiles(100);
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
        INFO("Using First 100 sequences of final.fasta and final.names") // Only appears on a FAIL
        
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
    
    SECTION("Remove Abundances") {
        INFO("Using GQY1XT001C44N8/size=3677/") // Only appears on a FAIL
        string seqName = "GQY1XT001C44N8/size=3677/";
        CAPTURE(testVParser.removeAbundances(seqName)); // Displays this variable on a FAIL
        
        CHECK(testVParser.removeAbundances(seqName) == "GQY1XT001C44N8");
    }
    
    
    SECTION("Create List File") {
        INFO("Using lines like: S	1	275	*	*	*	*	*	GQY1XT001C44N8/ab=3677/	*") // Only appears on a FAIL
        
        vsearchParser.getVsearchFile();
        ifstream in;
        testVParser.m->openInputFile(vsearchParser.getVsearchFile(), in);
        
        vector<string> seqNames;
        while (!in.eof()) {
            Sequence seq(in); testVParser.m->gobble(in);
            string name = seq.getName();
            seqNames.push_back(name);
        }
        in.close();
        testVParser.m->mothurRemove("tempSeqs.txt.sorted.fasta.temp");
        
        ofstream out;
        testVParser.m->openOutputFile("temp.txt", out);
        map<int, string> binNames;
        for (int i = 0; i < seqNames.size(); i++) {
            int bin = (i+1)%10;
            string name = testVParser.removeAbundances(seqNames[i]);
            //name = (testVParser.data.getNameMap())[name]; //dup names
            out << "S\t" + toString(bin) + "\t275\t*\t*\t*\t*\t*\t" + seqNames[i] + "\t*\n";
            
            map<int, string>::iterator it = binNames.find(bin);
            if (it != binNames.end()) { it->second += "," + name; }
            else { binNames[bin] = name; }
        }
        out.close();
        
        int numBins = binNames.size();
        
        testVParser.createListFile("temp.txt", "temp.list", "temp.rabund", "temp.sabund", numBins, "0.03");
        
        ifstream in2;
        testVParser.m->openInputFile("temp.list", in2);
        ListVector list(in2);
        in2.close();
        testVParser.m->mothurRemove("temp.list"); testVParser.m->mothurRemove("temp.rabund"); testVParser.m->mothurRemove("temp.sabund");
        
        //for each bin
        for (int i = 0; i < list.getNumBins(); i++) {
            string binnames = list.get(i);
            
            CAPTURE(binnames);
            
            CHECK(binnames == binNames[i]);
        }
    }
}
/**************************************************************************************************/
