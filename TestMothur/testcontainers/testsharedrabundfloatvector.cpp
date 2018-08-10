//
//  testsharedrabundfloatvector.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/9/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "testsharedrabundfloatvector.hpp"
#include "dataset.h"

/**************************************************************************************************/
TestSharedRabundFloatVector::TestSharedRabundFloatVector() : SharedRAbundFloatVector() {  //setup
    m = MothurOut::getInstance();
    
    TestDataSet data; relabundFile = data.getRelabundFile();
}
/**************************************************************************************************/
TestSharedRabundFloatVector::~TestSharedRabundFloatVector() {}//teardown
/**************************************************************************************************/

TEST(Test_Container_SharedRabundFloatVector, Constructors) {
    TestSharedRabundFloatVector test;
    
    SharedRAbundFloatVector noGroup;
    EXPECT_EQ(noGroup.getNumBins(), 0.0);
    EXPECT_EQ(noGroup.getNumSeqs(), 0.0);
    EXPECT_EQ(noGroup.getMaxRank(), 0.0);
    EXPECT_EQ(noGroup.getGroup(), "");
    
    SharedRAbundFloatVector noGroup2(10);
    EXPECT_EQ(noGroup2.getNumBins(), 10);
    EXPECT_EQ(noGroup2.getNumSeqs(), 0);
    EXPECT_EQ(noGroup2.getMaxRank(), 0);
    EXPECT_EQ(noGroup2.getGroup(), "");
    
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector noGroup3(abunds, 4.5, 10, 45.0);
    EXPECT_EQ(noGroup3.getNumBins(), 10);
    EXPECT_EQ(noGroup3.getNumSeqs(), 45.0);
    EXPECT_EQ(noGroup3.getMaxRank(), 4.5);
    EXPECT_EQ(noGroup3.getGroup(), "");
    
    SharedRAbundFloatVector noGroup4(abunds);
    EXPECT_EQ(noGroup3.getNumBins(), 10);
    EXPECT_EQ(noGroup3.getNumSeqs(), 45.0);
    EXPECT_EQ(noGroup3.getMaxRank(), 4.5);
    EXPECT_EQ(noGroup3.getGroup(), "");
    
    SharedRAbundFloatVector noGroup5(noGroup4);
    EXPECT_EQ(noGroup3.getNumBins(), 10);
    EXPECT_EQ(noGroup3.getNumSeqs(), 45.0);
    EXPECT_EQ(noGroup3.getMaxRank(), 4.5);
    EXPECT_EQ(noGroup3.getGroup(), "");
    
    ifstream in;
    Utils util; util.openInputFile(test.relabundFile, in);
    util.getline(in); //gobble headers
    SharedRAbundFloatVector temp(in);
    EXPECT_EQ(temp.getNumBins(), 58);
    ASSERT_NEAR(temp.getNumSeqs(), 1, 0.001);
    ASSERT_NEAR(temp.getMaxRank(), 0.25, 0.001);
    EXPECT_EQ(temp.getGroup(), "F003D000");
    
    SharedRAbundFloatVector temp2(in);
    EXPECT_EQ(temp2.getNumBins(), 58);
    ASSERT_NEAR(temp2.getNumSeqs(), 1, 0.001);
    ASSERT_NEAR(temp2.getMaxRank(), 0.16667, 0.001);
    EXPECT_EQ(temp2.getGroup(), "F003D002");
    
    int numBins; string label, groupN;
    in >> label >> groupN >> numBins;
    SharedRAbundFloatVector temp3(in, label, groupN, numBins);
    EXPECT_EQ(temp3.getNumBins(), 58);
    ASSERT_NEAR(temp3.getNumSeqs(), 1, 0.001);
    ASSERT_NEAR(temp3.getMaxRank(), 0.277778, 0.001);
    EXPECT_EQ(temp3.getGroup(), "F003D004");
    in.close();
}

TEST(Test_Container_SharedRabundFloatVector, GetsSets) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getGroup(), "");
    
    temp.set(5, 0.5);
    EXPECT_EQ(temp.getNumSeqs(), 41.5);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.get()[5], 0.5);
    EXPECT_EQ(temp.get(5), 0.5);
    EXPECT_EQ(temp.getGroup(), "");
    temp.setGroup("myGroup");
    EXPECT_EQ(temp.getGroup(), "myGroup");
}

TEST(Test_Container_SharedRabundFloatVector, PushBack) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getGroup(), "");

    temp.push_back(0.25);
    EXPECT_EQ(temp.getNumSeqs(), 45.25);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getNumBins(), 11);
}

TEST(Test_Container_SharedRabundFloatVector, ResizeSize) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getGroup(), "");
    
    temp.resize(20);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getNumBins(), 20);
    
    temp.resize(5);
    EXPECT_EQ(temp.getNumSeqs(), 22.5);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getNumBins(), 5);
    EXPECT_EQ(temp.size(), 5);
}

TEST(Test_Container_SharedRabundFloatVector, ClearRemove) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getGroup(), "");
    
    
    temp.remove(2);
    EXPECT_EQ(temp.getNumSeqs(), 40.5);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getNumBins(), 9);
    
    temp.clear();
    EXPECT_EQ(temp.getNumSeqs(), 0);
    EXPECT_EQ(temp.getMaxRank(), 0);
    EXPECT_EQ(temp.getNumBins(), 0);
    EXPECT_EQ(temp.size(), 0);
}

TEST(Test_Container_SharedRabundFloatVector, GetRabundVector) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    EXPECT_EQ(temp.getGroup(), "");
    
    RAbundVector rabund = temp.getRAbundVector();
    EXPECT_EQ(rabund.getNumBins(), 10);
    EXPECT_EQ(rabund.getNumSeqs(), 40);
    EXPECT_EQ(rabund.getMaxRank(), 4);
    EXPECT_EQ(rabund.get(5), 4);
    
}

TEST(Test_Container_SharedRabundFloatVector, GetSabundVector) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    
    SAbundVector sabund = temp.getSAbundVector();
    EXPECT_EQ(sabund.getNumBins(), 10);
    EXPECT_EQ(sabund.getNumSeqs(), 40);
    EXPECT_EQ(sabund.getMaxRank(), 4);
    EXPECT_EQ(sabund.get(4), 10);
    
}

TEST(Test_Container_SharedRabundFloatVector, RAbundFloatVector) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    
    RAbundFloatVector rabundFloat = temp.getRAbundFloatVector();
    EXPECT_EQ(rabundFloat.getNumBins(), 10.0);
    EXPECT_EQ(rabundFloat.getNumSeqs(), 45.0);
    EXPECT_EQ(rabundFloat.getMaxRank(), 4.5);
    rabundFloat.set(5, 12.5);
    EXPECT_EQ(rabundFloat.getNumBins(), 10.0);
    EXPECT_EQ(rabundFloat.getNumSeqs(), 53.0);
    EXPECT_EQ(rabundFloat.getMaxRank(), 12.5);
    
}

TEST(Test_Container_SharedRabundFloatVector, SharedRAbundVector) {
    vector<float> abunds(10, 4.5);
    SharedRAbundFloatVector temp(abunds);
    EXPECT_EQ(temp.getNumBins(), 10);
    EXPECT_EQ(temp.getNumSeqs(), 45.0);
    EXPECT_EQ(temp.getMaxRank(), 4.5);
    
    SharedRAbundVector rabund = temp.getSharedRAbundVector();
    EXPECT_EQ(rabund.getNumBins(), 10);
    EXPECT_EQ(rabund.getNumSeqs(), 40);
    EXPECT_EQ(rabund.getMaxRank(), 4);
    rabund.set(5, 12);
    EXPECT_EQ(rabund.getNumBins(), 10);
    EXPECT_EQ(rabund.getNumSeqs(), 48);
    EXPECT_EQ(rabund.getMaxRank(), 12);
}

/**************************************************************************************************/

