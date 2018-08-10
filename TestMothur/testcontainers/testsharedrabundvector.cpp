//
//  testsharedrabundvector.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/9/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testsharedrabundvector.hpp"
#include "dataset.h"

/**************************************************************************************************/
TestSharedRabundVector::TestSharedRabundVector() : SharedRAbundVector() {  //setup
    m = MothurOut::getInstance();
    
    TestDataSet data; sharedFile = data.getSharedFile();
}
/**************************************************************************************************/
TestSharedRabundVector::~TestSharedRabundVector() {}//teardown
/**************************************************************************************************/

TEST(Test_Container_SharedRabundVector, Constructors) {
    TestSharedRabundVector test;
    
    SharedRAbundVector noGroup;
    EXPECT_EQ(noGroup.getNumBins(), 0);
    EXPECT_EQ(noGroup.getNumSeqs(), 0);
    EXPECT_EQ(noGroup.getMaxRank(), 0);
    EXPECT_EQ(noGroup.getGroup(), "");
    
    SharedRAbundVector noGroup2(10);
    EXPECT_EQ(noGroup2.getNumBins(), 10);
    EXPECT_EQ(noGroup2.getNumSeqs(), 0);
    EXPECT_EQ(noGroup2.getMaxRank(), 0);
    EXPECT_EQ(noGroup2.getGroup(), "");
    
    vector<int> abunds(10, 5);
    SharedRAbundVector noGroup3(abunds, 5, 10, 50);
    EXPECT_EQ(noGroup3.getNumBins(), 10);
    EXPECT_EQ(noGroup3.getNumSeqs(), 50);
    EXPECT_EQ(noGroup3.getMaxRank(), 5);
    EXPECT_EQ(noGroup3.getGroup(), "");
    
    SharedRAbundVector noGroup4(abunds);
    EXPECT_EQ(noGroup4.getNumBins(), 10);
    EXPECT_EQ(noGroup4.getNumSeqs(), 50);
    EXPECT_EQ(noGroup4.getMaxRank(), 5);
    EXPECT_EQ(noGroup4.getGroup(), "");
    
    SharedRAbundVector noGroup5(noGroup4);
    EXPECT_EQ(noGroup5.getNumBins(), 10);
    EXPECT_EQ(noGroup5.getNumSeqs(), 50);
    EXPECT_EQ(noGroup5.getMaxRank(), 5);
    EXPECT_EQ(noGroup5.getGroup(), "");
    
    ifstream in;
    Utils util; util.openInputFile(test.sharedFile, in);
    util.getline(in); //gobble headers
    SharedRAbundVector temp(in);
    EXPECT_EQ(temp.getNumBins(), 58);
    EXPECT_EQ(temp.getNumSeqs(), 20);
    EXPECT_EQ(temp.getMaxRank(), 5);
    EXPECT_EQ(temp.getGroup(), "F003D000");
    
    SharedRAbundVector temp2(in);
    EXPECT_EQ(temp2.getNumBins(), 58);
    EXPECT_EQ(temp2.getNumSeqs(), 18);
    EXPECT_EQ(temp2.getMaxRank(), 3);
    EXPECT_EQ(temp2.getGroup(), "F003D002");
    
    int numBins; string label, groupN;
    in >> label >> groupN >> numBins;
    SharedRAbundVector temp3(in, label, groupN, numBins);
    EXPECT_EQ(temp3.getNumBins(), 58);
    EXPECT_EQ(temp3.getNumSeqs(), 18);
    EXPECT_EQ(temp3.getMaxRank(), 5);
    EXPECT_EQ(temp3.getGroup(), "F003D004");
    in.close();
}

TEST(Test_Container_SharedRabundVector, GetsSets) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    EXPECT_EQ(temp.getGroup(), "");
    
    temp.set(5, 20);
    EXPECT_EQ(temp.getNumSeqs(), 65);
    EXPECT_EQ(temp.getMaxRank(), 20);
    EXPECT_EQ(temp.get()[5], 20);
    EXPECT_EQ(temp.get(5), 20);
    EXPECT_EQ(temp.getGroup(), "");
    temp.setGroup("myGroup");
    EXPECT_EQ(temp.getGroup(), "myGroup");
}

TEST(Test_Container_SharedRabundVector, Increment) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    
    temp.increment(5);
    EXPECT_EQ(temp.getNumSeqs(), 51);
    EXPECT_EQ(temp.getMaxRank(), 6);
}

TEST(Test_Container_SharedRabundVector, PushBack) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    
    temp.push_back(20);
    EXPECT_EQ(temp.getNumSeqs(), 70);
    EXPECT_EQ(temp.getMaxRank(), 20);
    EXPECT_EQ(temp.getNumBins(), 11);
}

TEST(Test_Container_SharedRabundVector, ResizeSize) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    
    temp.resize(20);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    EXPECT_EQ(temp.getNumBins(), 20);
    
    temp.resize(5);
    EXPECT_EQ(temp.getNumSeqs(), 25);
    EXPECT_EQ(temp.getMaxRank(), 5);
    EXPECT_EQ(temp.getNumBins(), 5);
    EXPECT_EQ(temp.size(), 5);
}

TEST(Test_Container_SharedRabundVector, ClearRemove) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    
    temp.remove(2);
    EXPECT_EQ(temp.getNumSeqs(), 45);
    EXPECT_EQ(temp.getMaxRank(), 5);
    EXPECT_EQ(temp.getNumBins(), 9);
    
    temp.clear();
    EXPECT_EQ(temp.getNumSeqs(), 0);
    EXPECT_EQ(temp.getMaxRank(), 0);
    EXPECT_EQ(temp.getNumBins(), 0);
    EXPECT_EQ(temp.size(), 0);
}

TEST(Test_Container_SharedRabundVector, GetRabundVector) {
    vector<int> abunds(10, 10);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 100);
    EXPECT_EQ(temp.getMaxRank(), 10);
    
    RAbundVector rabund = temp.getRAbundVector();
    EXPECT_EQ(rabund.getNumBins(), 10);
    EXPECT_EQ(rabund.getNumSeqs(), 100);
    EXPECT_EQ(rabund.getMaxRank(), 10);
}

TEST(Test_Container_SharedRabundVector, GetSabundVector) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    
    SAbundVector sabund = temp.getSAbundVector();
    EXPECT_EQ(sabund.getNumBins(), 10);
    EXPECT_EQ(sabund.getNumSeqs(), 50);
    EXPECT_EQ(sabund.getMaxRank(), 5);
    EXPECT_EQ(sabund.get(5), 10);
}

TEST(Test_Container_SharedRabundVector, RAbundFloatVector) {
    vector<int> abunds(10, 5);
    SharedRAbundVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 50);
    EXPECT_EQ(temp.getMaxRank(), 5);
    
    RAbundFloatVector rabundFloat = temp.getRAbundFloatVector();
    EXPECT_EQ(rabundFloat.getNumBins(), 10.0);
    EXPECT_EQ(rabundFloat.getNumSeqs(), 50.0);
    EXPECT_EQ(rabundFloat.getMaxRank(), 5.0);
    rabundFloat.set(5, 12.5);
    EXPECT_EQ(rabundFloat.getNumBins(), 10.0);
    EXPECT_EQ(rabundFloat.getNumSeqs(), 57.5);
    EXPECT_EQ(rabundFloat.getMaxRank(), 12.5);
}

/**************************************************************************************************/



