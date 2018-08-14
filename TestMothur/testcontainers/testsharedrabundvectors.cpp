//
//  testsharedrabundvectors.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testsharedrabundvectors.hpp"
#include "dataset.h"

/**************************************************************************************************/
TestSharedRabundVectors::TestSharedRabundVectors() : SharedRAbundVectors() {  //setup
    m = MothurOut::getInstance();
    
    TestDataSet data; sharedFile = data.getSharedFile();
}
/**************************************************************************************************/
TestSharedRabundVectors::~TestSharedRabundVectors() {}//teardown
/**************************************************************************************************/

TEST(Test_Container_SharedRabundVectors, Constructors) {
    TestSharedRabundVectors test;
    
    SharedRAbundVectors temp;
    EXPECT_EQ(temp.getNumBins(), 0);
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    EXPECT_EQ(nextLabel, "");
    EXPECT_EQ(labelTag, "Otu");
    EXPECT_EQ(userGroups[0], "F003D000");
    
    SharedRAbundVectors copy(fileRead);
    EXPECT_EQ(copy.getLabel(), "0.03");
    EXPECT_EQ(labelTag, "Otu");
    EXPECT_EQ(userGroups[0], "F003D000");

}

TEST(Test_Container_SharedRabundVectors, GetsSets) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    EXPECT_EQ(nextLabel, "");
    EXPECT_EQ(labelTag, "Otu");
    EXPECT_EQ(userGroups[0], "F003D000");
    
    //setLabels
    fileRead.setLabels("0.99");
    EXPECT_EQ(fileRead.getLabel(), "0.99");
    
    //getOTUTotal
    EXPECT_EQ(fileRead.getOTUTotal(4), 12);
    
    //getOTU
    vector<int> otu5 = fileRead.getOTU(4);
    EXPECT_EQ(otu5[0], 0);
    EXPECT_EQ(otu5[1], 2);
    
    //get
    EXPECT_EQ(fileRead.get(4, "F003D000"), 0);
    EXPECT_EQ(fileRead.get(4, "F003D002"), 2);
    EXPECT_EQ(fileRead.get(0, "F003D142"), 6);
    
    //set
    fileRead.set(4, 10, "F003D000");
    EXPECT_EQ(fileRead.get(4, "F003D000"), 10);
    fileRead.set(4, 15, "F003D002");
    EXPECT_EQ(fileRead.get(4, "F003D002"), 15);
    fileRead.set(0, 3, "F003D142");
    EXPECT_EQ(fileRead.get(0, "F003D142"), 3);
    
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
    EXPECT_EQ(fileRead.getNumSeqsSmallestGroup(), 15);
    
    //getNamesGroups
    EXPECT_EQ(fileRead.getNamesGroups()[1], "F003D002");
    
    //getNumGroups
    EXPECT_EQ(fileRead.getNumGroups(), 10);
    
    //getNumSeqs
    EXPECT_EQ(fileRead.getNumSeqs("F003D002"), 31);
    
}

TEST(Test_Container_SharedRabundVectors, PushBack) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    //getNumGroups
    EXPECT_EQ(fileRead.getNumGroups(), 10);
    
    vector<int> abunds(58, 5);
    SharedRAbundVector* temp = new SharedRAbundVector(abunds);
    temp->setGroup("myabunds");
    fileRead.push_back(temp);
    
    //getNumGroups
    EXPECT_EQ(fileRead.getNumGroups(), 11);
    
    vector<int> otuAbunds(11, 2); string otuLabel = "Otu59";
    fileRead.push_back(otuAbunds, otuLabel);
    EXPECT_EQ(fileRead.getNumBins(), 59);
}

TEST(Test_Container_SharedRabundVectors, eliminateZeroOTUS) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    EXPECT_EQ(fileRead.getNumBins(), 58);
    
    fileRead.set(16, 0, "F003D142");
    EXPECT_EQ(fileRead.get(16, "F003D142"), 0); //zero out bin
    fileRead.eliminateZeroOTUS(); //remove zeroed out bin
    
    EXPECT_EQ(fileRead.getNumBins(), 57);
    
}

TEST(Test_Container_SharedRabundVectors, Removes) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    EXPECT_EQ(fileRead.getNumBins(), 58);
    EXPECT_EQ(fileRead.removeOTU(16), 2);
    EXPECT_EQ(fileRead.getNumBins(), 57);
    
    vector<string> groups;
    groups.push_back("F003D142");
    groups.push_back("F003D002");
    groups.push_back("F003D004");
    groups.push_back("F003D006");
    fileRead.removeGroups(groups);
    
    EXPECT_EQ(fileRead.getNumBins(), 41);
    EXPECT_EQ(fileRead.getNumSeqsSmallestGroup(), 15);
    fileRead.removeGroups(20); //remove groups with abundance less than 20
    EXPECT_EQ(fileRead.getNumGroups(), 4);
}

TEST(Test_Container_SharedRabundVectors, SizeClear) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    EXPECT_EQ(fileRead.size(), 10);
    fileRead.clear();
    EXPECT_EQ(fileRead.getNumBins(), 0);
}

TEST(Test_Container_SharedRabundVectors, GetRabundVector) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    
    RAbundVector temp = fileRead.getRAbundVector();
    EXPECT_EQ(temp.get(0), 24);
    
    RAbundVector temp2 = fileRead.getRAbundVector("F003D142");
    EXPECT_EQ(temp2.get(0), 6);
}

TEST(Test_Container_SharedRabundVectors, GetSabundVector) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    
    SAbundVector temp = fileRead.getSAbundVector();
    EXPECT_EQ(temp.get(5), 1); //number of OTUs with abundance of 5
    EXPECT_EQ(temp.get(1), 35); //number of OTUs with abundance of 1
    EXPECT_EQ(temp.getMaxRank(), 24);
}

TEST(Test_Container_SharedRabundVectors, GetSharedRAbundVectors) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    
    vector<SharedRAbundVector*> temp = fileRead.getSharedRAbundVectors();
    EXPECT_EQ(temp[0]->get(5), 0); //first groups abundance of OTU5
    EXPECT_EQ(temp[1]->get(5), 2); //first groups abundance of OTU5
    EXPECT_EQ(temp[0]->get(0), 1); //first groups abundance of OTU1
    
}

TEST(Test_Container_SharedRabundVectors, GetSharedRAbundFloatVectors) {
    TestSharedRabundVectors test;
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    vector<string> userGroups; string nextLabel, labelTag;
    SharedRAbundVectors fileRead(in, userGroups, nextLabel, labelTag);
    
    
    vector<SharedRAbundFloatVector*> temp = fileRead.getSharedRAbundFloatVectors();
    EXPECT_EQ(temp[0]->get(5), 0.0); //first groups abundance of OTU5 as a float. This is not the same as if we ran get.relabund() command
    EXPECT_EQ(temp[1]->get(5), 2.0); //first groups abundance of OTU5 as a float. This is not the same as if we ran get.relabund() command
    EXPECT_EQ(temp[0]->get(0), 1.0); //first groups abundance of OTU1 as a float. This is not the same as if we ran get.relabund() command
    
}


/**************************************************************************************************/


