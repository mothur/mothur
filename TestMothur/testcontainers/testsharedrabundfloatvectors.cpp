//
//  testsharedrabundfloatvectors.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/14/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testsharedrabundfloatvectors.hpp"
#include "dataset.h"

/**************************************************************************************************/
TestSharedRabundFloatVectors::TestSharedRabundFloatVectors() : SharedRAbundFloatVectors() {  //setup
    m = MothurOut::getInstance();
    
    TestDataSet data; relabundFile = data.getRelabundFile();
}
/**************************************************************************************************/
TestSharedRabundFloatVectors::~TestSharedRabundFloatVectors() {}//teardown
/**************************************************************************************************/

TEST(Test_Container_SharedRabundFloatVectors, Constructors) {
    TestSharedRabundFloatVectors test;
    
    SharedRAbundFloatVectors temp;
    EXPECT_EQ(temp.getNumBins(), 0);
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    EXPECT_EQ(nextLabel, "");
    EXPECT_EQ(labelTag, "Otu");
    EXPECT_EQ(userGroups[0], "F003D000");
    
    SharedRAbundFloatVectors copy(fileRead);
    EXPECT_EQ(copy.getLabel(), "0.03");
    EXPECT_EQ(labelTag, "Otu");
    EXPECT_EQ(userGroups[0], "F003D000");
    
}

TEST(Test_Container_SharedRabundFloatVectors, GetsSets) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    EXPECT_EQ(nextLabel, "");
    EXPECT_EQ(labelTag, "Otu");
    EXPECT_EQ(userGroups[0], "F003D000");
    
    //setLabels
    fileRead.setLabels("0.99");
    EXPECT_EQ(fileRead.getLabel(), "0.99");
    
    //getOTUTotal
    ASSERT_NEAR(fileRead.getOTUTotal(4), 0.615171, 0.001);
    
    //getOTU
    vector<float> otu5 = fileRead.getOTU(4);
    ASSERT_NEAR(otu5[0], 0.0, 0.001);
    ASSERT_NEAR(otu5[2], 0.055556, 0.001);
    
    //get
    ASSERT_NEAR(fileRead.get(4, "F003D000"), 0.0, 0.001);
    ASSERT_NEAR(fileRead.get(4, "F003D002"), 0.111111, 0.001);
    ASSERT_NEAR(fileRead.get(0, "F003D142"), 0.222222, 0.001);
    
    //set
    fileRead.set(4, 0.25, "F003D000");
    ASSERT_NEAR(fileRead.get(4, "F003D000"), 0.25, 0.001);
    fileRead.set(4, 0.1234, "F003D002");
    ASSERT_NEAR(fileRead.get(4, "F003D002"), 0.1234, 0.001);
    fileRead.set(0, 0.456, "F003D142");
    ASSERT_NEAR(fileRead.get(0, "F003D142"), 0.456, 0.001);
    
    //getOTUNames
    vector<string> otuNames = fileRead.getOTUNames();
    EXPECT_EQ(otuNames[10], "Otu11");
    
    //setOTUNames
    otuNames[5] = "Otu99";
    fileRead.setOTUNames(otuNames);
    EXPECT_EQ(fileRead.getOTUNames()[5], "Otu99");
    EXPECT_EQ(fileRead.getOTUName(5), "Otu99");
    fileRead.setOTUName(5, "Otu6");
    EXPECT_EQ(fileRead.getOTUName(5), "Otu6");
    
    //getNumBins
    EXPECT_EQ(fileRead.getNumBins(), 58);
    
    //getNumSeqsSmallestGroup
    ASSERT_NEAR(fileRead.getNumSeqsSmallestGroup(), 1.0, 0.001);
    
    //getNamesGroups
    EXPECT_EQ(fileRead.getNamesGroups()[1], "F003D002");
    
    //getNumGroups
    EXPECT_EQ(fileRead.getNumGroups(), 10);
    
    //getNumSeqs
    ASSERT_NEAR(fileRead.getNumSeqs("F003D002"), 1.1234, 0.001); //1.0 + 0.1234 (set from above)
}

TEST(Test_Container_SharedRabundFloatVectors, PushBack) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    //getNumGroups
    EXPECT_EQ(fileRead.getNumGroups(), 10);
    
    vector<float> abunds(58, 0.5);
    SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(abunds);
    temp->setGroup("myabunds");
    fileRead.push_back(temp);
    
    //getNumGroups
    EXPECT_EQ(fileRead.getNumGroups(), 11);
}


TEST(Test_Container_SharedRabundFloatVectors, eliminateZeroOTUS) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    EXPECT_EQ(fileRead.getNumBins(), 58);
    
    fileRead.set(16, 0, "F003D142");
    EXPECT_EQ(fileRead.get(16, "F003D142"), 0); //zero out bin
    fileRead.eliminateZeroOTUS(); //remove zeroed out bin
    
    EXPECT_EQ(fileRead.getNumBins(), 57);
    
}

TEST(Test_Container_SharedRabundFloatVectors, Removes) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    EXPECT_EQ(fileRead.getNumBins(), 58);
    ASSERT_NEAR(fileRead.removeOTU(16), 0.074074, 0.001);
    EXPECT_EQ(fileRead.getNumBins(), 57);
    
    vector<string> groups;
    groups.push_back("F003D142");
    groups.push_back("F003D002");
    groups.push_back("F003D004");
    groups.push_back("F003D006");
    fileRead.removeGroups(groups);
    
    EXPECT_EQ(fileRead.getNumBins(), 41);
    ASSERT_NEAR(fileRead.getNumSeqsSmallestGroup(), 0.999999, 0.001);
    fileRead.removeGroups(20); //remove groups with abundance less than 20
    EXPECT_EQ(fileRead.getNumGroups(), 0);
}

TEST(Test_Container_SharedRabundFloatVectors, SizeClear) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    EXPECT_EQ(fileRead.size(), 10);
    fileRead.clear();
    EXPECT_EQ(fileRead.getNumBins(), 0);
}

TEST(Test_Container_SharedRabundFloatVectors, GetRabundVector) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    
    RAbundVector temp = fileRead.getRAbundVector();
    EXPECT_EQ(temp.get(0), 1);
    
}

TEST(Test_Container_SharedRabundFloatVectors, GetSabundVector) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    SAbundVector temp = fileRead.getSAbundVector();
    EXPECT_EQ(temp.get(0), 0); //number of OTUs with abundance of 5
    EXPECT_EQ(temp.get(1), 1); //number of OTUs with abundance of 1
    EXPECT_EQ(temp.getMaxRank(), 1);
}

TEST(Test_Container_SharedRabundFloatVectors, GetSharedRAbundVectors) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    vector<SharedRAbundVector*> temp = fileRead.getSharedRAbundVectors();
    EXPECT_EQ(temp[0]->get(5), 0); //first groups abundance of OTU5
    EXPECT_EQ(temp[1]->get(5), 0); //first groups abundance of OTU5
    EXPECT_EQ(temp[0]->get(0), 0); //first groups abundance of OTU1
}

TEST(Test_Container_SharedRabundFloatVectors, GetSharedRAbundFloatVectors) {
    TestSharedRabundFloatVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundFloatVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    vector<SharedRAbundFloatVector*> temp = fileRead.getSharedRAbundFloatVectors();
    EXPECT_EQ(temp[0]->get(5), 0.0); //first groups abundance of OTU5 as a float.
    ASSERT_NEAR(temp[1]->get(5), 0.111111, 0.001); //first groups abundance of OTU5 as a float.
    ASSERT_NEAR(temp[0]->get(0), 0.05, 0.001); //first groups abundance of OTU1 as a float.
}
/**************************************************************************************************/


