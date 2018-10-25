//
//  testcounttable.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/25/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testcounttable.hpp"

/**************************************************************************************************/
TestCountTable::TestCountTable() {  //setup
    m = MothurOut::getInstance();
    
    TestDataSet data;
    vector<string> filenames = data.getSubsetFNGFiles();
    fastafile = filenames[0];
    namefile = filenames[1];
    groupfile = filenames[2];
    countfile = data.getCountTableFile();
}
/**************************************************************************************************/
TestCountTable::~TestCountTable() {}//teardown
/**************************************************************************************************/
//Testing createTable functions
TEST(Test_Container_CountTable, createTables) {
    //int createTable(string, string, bool); //namefile, groupfile, createGroup
    TestCountTable testData;
    CountTable ct;
    ct.createTable(testData.namefile, testData.groupfile, false);
    
    EXPECT_EQ(ct.getNumGroups(), 10);
    EXPECT_EQ(ct.getNumSeqs(), 200);
    
    //int createTable(set<string>&, map<string, string>&, set<string>&); //seqNames, seqName->group, groupNames
    set<string> seqNames;
    seqNames.insert("seq1"); seqNames.insert("seq2"); seqNames.insert("seq3"); seqNames.insert("seq4"); seqNames.insert("seq5");
    set<string> groupNames;
    groupNames.insert("group1"); groupNames.insert("group2");
    map<string, string> groupMap;
    groupMap["seq1"] = "group1"; groupMap["seq2"] = "group1"; groupMap["seq3"] = "group1"; groupMap["seq4"] = "group2"; groupMap["seq5"] = "group2";
    
    ct.clearTable();
    ct.createTable(seqNames, groupMap, groupNames);
    
    EXPECT_EQ(ct.getNumGroups(), 2);
    EXPECT_EQ(ct.getNumSeqs(), 5);
    
    //int readTable(string, bool, bool); //filename, readGroups, mothurRunning
    ct.clearTable();
    ct.readTable(testData.countfile, true, true);
    
    EXPECT_EQ(ct.getNumGroups(), 10);
    EXPECT_EQ(ct.getNumSeqs(), 200);
    
    ct.clearTable();
    ct.readTable(testData.countfile, false, true);
    
    EXPECT_EQ(ct.getNumGroups(), 0);
    EXPECT_EQ(ct.getNumSeqs(), 200);

    //int readTable(string, string); //filename, format - if format=fasta, read fasta file and create unique table
    ct.clearTable();
    ct.readTable(testData.fastafile, "fasta");
    
    EXPECT_EQ(ct.getNumGroups(), 0);
    EXPECT_EQ(ct.getNumSeqs(), 93);
    
    EXPECT_EQ(ct.getHardCodedHeaders()[0], "Representative_Sequence");
    EXPECT_EQ(ct.getHardCodedHeaders()[1], "total");
}
/**************************************************************************************************/
//Testing testGroups functions
TEST(Test_Container_CountTable, testGroups) {
    TestCountTable testData;
    CountTable ct;
    EXPECT_EQ(ct.testGroups(testData.countfile), true);
    
    vector<string> groups;
    ct.testGroups(testData.countfile, groups);
    
    EXPECT_EQ(groups[0], "F003D000");
    EXPECT_EQ(groups[1], "F003D002");
    EXPECT_EQ(groups[2], "F003D004");
    EXPECT_EQ(groups[3], "F003D006");
    EXPECT_EQ(groups[4], "F003D008");
    EXPECT_EQ(groups[5], "F003D142");
    
    ct.createTable(testData.namefile, testData.groupfile, false);
    CountTable ct2; ct2.copy(&ct);
    EXPECT_EQ(ct2.getNumGroups(), 10);
    EXPECT_EQ(ct2.getNumSeqs(), 200);
    EXPECT_EQ(ct2.hasGroupInfo(), true);
    
    groups = ct2.getNamesOfGroups();
    EXPECT_EQ(groups[0], "F003D000");
    EXPECT_EQ(groups[1], "F003D002");
    EXPECT_EQ(groups[2], "F003D004");
    EXPECT_EQ(groups[3], "F003D006");
    EXPECT_EQ(groups[4], "F003D008");
    EXPECT_EQ(groups[5], "F003D142");
    
    groups.clear();
    groups.push_back("group1"); groups.push_back("group2"); groups.push_back("group3");
    ct2.setNamesOfGroups(groups);
    groups = ct2.getNamesOfGroups();
    EXPECT_EQ(groups[0], "group1");
    EXPECT_EQ(groups[1], "group2");
    EXPECT_EQ(groups[2], "group3");
    
    ct2.addGroup("group4");
    ct2.removeGroup("group2");
    groups = ct2.getNamesOfGroups();
    EXPECT_EQ(groups[0], "group1");
    EXPECT_EQ(groups[1], "group3");
    EXPECT_EQ(groups[2], "group4");
}
/**************************************************************************************************/
//Testing testGroups functions
TEST(Test_Container_CountTable, push_backs) {
    TestCountTable testData;
    CountTable ct;
    ct.push_back("seq1");
    
    EXPECT_EQ(ct.getNamesOfSeqs()[0], "seq1");
    
    ct.push_back("seq2", 15);
    EXPECT_EQ(ct.getNumSeqs(), 16);
    EXPECT_EQ(ct.size(), 2);
    
    ct.renameSeq("seq1", "mySeq");
    EXPECT_EQ(ct.getNamesOfSeqs()[0], "mySeq");
    
    ct.remove("mySeq");
    EXPECT_EQ(ct.getNamesOfSeqs()[0], "seq2");
    ct.push_back("seq3", 10);
    EXPECT_EQ(ct.get("seq2"), 0);
    EXPECT_EQ(ct.get("seq3"), 1);
}
/**************************************************************************************************/
//Testing testGroups functions
TEST(Test_Container_CountTable, push_backGroups) {
    TestCountTable testData;
    CountTable ct;
    ct.createTable(testData.namefile, testData.groupfile, false);
    ct.setAbund("GQY1XT001B1CEF", "F003D000", 50);
    EXPECT_EQ(ct.getGroupCount("GQY1XT001B1CEF", "F003D000"), 50);
    
    vector<int> abunds; abunds.resize(10, 100);
    ct.push_back("mySeq", abunds);
    EXPECT_EQ(ct.getGroupCount("mySeq", "F003D000"), 100);
    
}
/**************************************************************************************************/

//vector<string> getGroups(string); //returns vector of groups represented by this sequences
//vector<int> getGroupCounts(string);  //returns group counts for a seq passed in, if no group info is in file vector is blank. Order is the same as the groups returned by getGroups function.
//int getGroupCount(string, string); //returns number of seqs for that group for that seq
//int getGroupCount(string); // returns total seqs for that group
//int getNumSeqs(string); //returns total seqs for that seq, 0 if not found
//int setNumSeqs(string, int); //set total seqs for that seq, return -1 if not found
//int getNumSeqs() { return total; } //return total number of seqs
//int getNumUniqueSeqs() { return uniques; } //return number of unique/representative seqs
//int getGroupIndex(string); //returns index in getGroupCounts vector of specific group

//vector<string> getNamesOfSeqs(); //return names of all seqeunce in table
//vector<string> getNamesOfSeqs(string); //returns names of seqs in specific group in table
//int mergeCounts(string, string); //combines counts for 2 seqs, saving under the first name passed in.
//ListVector getListVector();
//SharedRAbundVectors* getShared();
//SharedRAbundVectors* getShared(vector<string>); //set of groups selected
//map<string, int> getNameMap();  //sequenceName -> total number of sequences it represents



