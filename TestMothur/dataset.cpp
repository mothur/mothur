//
//  dataset.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/24/16.
//  Copyright (c) 2016 Schloss Lab. All rights reserved.
//

#include "dataset.h"
#include "inputdata.h"

/***********************************************************************/
TestDataSet::TestDataSet() {
    m = MothurOut::getInstance();
    gMap = NULL;
    testDir = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/";
}
/***********************************************************************/
void TestDataSet::createCountTable() {
    fillGroup();
    fillNames();
    ct = new CountTable();
    
    for (map<string, string>::iterator itNameMap = nameMap.begin(); itNameMap !=nameMap.end(); itNameMap++) {
        string firstCol = itNameMap->first;
        string secondCol = itNameMap->second;
        vector<string> names;
        util.splitAtChar(secondCol, names, ',');
        
        //set to 0
        map<string, int> groupCounts;
        int total = 0;
        vector<string> Groups = gMap->getNamesOfGroups();
        ct->setNamesOfGroups(Groups);
        for (int i = 0; i < Groups.size(); i++) { groupCounts[Groups[i]] = 0; }
        
        //get counts for each of the users groups
        for (int i = 0; i < names.size(); i++) {
            string group = gMap->getGroup(names[i]);
            
            map<string, int>::iterator it = groupCounts.find(group);
            
            //if not found, then this sequence is not from a group we care about
            if (it != groupCounts.end()) {
                it->second++;
                total++;
            }
        }
        
        if (total != 0) {
            vector<int> abunds;
            for (map<string, int>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) { abunds.push_back(it->second); }
            ct->push_back(firstCol, abunds);
        }
    }
    delete gMap; gMap = NULL;
    nameMap.clear();
}
/***********************************************************************/

vector<string> TestDataSet::getSubsetFNGFiles() {
    vector<string> filenames;
    filenames.push_back(testDir+"test.fasta");
    filenames.push_back(testDir+"test.names");
    filenames.push_back(testDir+"test.groups");
    
    return filenames;
}
/***********************************************************************/

vector<string> TestDataSet::getOptiRefFiles() {
    vector<string> filenames;
    filenames.push_back(testDir+"silva.v4.unique.fasta");
    filenames.push_back(testDir+"silva.v4.count_table");
    filenames.push_back(testDir+"silva.v4.unique.dist");
    filenames.push_back(testDir+"silva.v4.unique.phylip.dist");
    filenames.push_back(testDir+"silva.v4.unique.opti_mcc.list");
    filenames.push_back(testDir+"test.fit.dist");
    
    return filenames;
}

/***********************************************************************/

vector<string> TestDataSet::getOligosFiles() {
    vector<string> filenames;
    filenames.push_back(testDir+"GQY1XT001.oligos"); //single
    filenames.push_back(testDir+"Undetermined.oligos"); //paired
    filenames.push_back(testDir+"bo.oligos"); //index "NONE"
    filenames.push_back(testDir+"comboNames.oligos"); //named primers and named barcodes
    
    return filenames;
}
/***********************************************************************/
string TestDataSet::getSubsetFNGDistFile() {
    
    return (testDir+"test.dist");
}
/***********************************************************************/
string TestDataSet::getSubsetFNGPhylipDistFile() {
    
    return (testDir+"test.phylip.dist");
}
/***********************************************************************/
void TestDataSet::fillSeqs() {
    seqs.clear();
    
    //read info from stable file
    //string testfile = m->getTestFilePath() + "testFile.fasta";
     string testfile = testDir+"test.fasta";
    
    ifstream in;
    util.openInputFile(testfile, in);
    
    while (!in.eof()) {
        if (m->getControl_pressed()) { break; }
        
        Sequence read(in); util.gobble(in);
        seqs.push_back(read);
    }
    in.close();
}
/***********************************************************************/
void TestDataSet::fillNames() {
    nameMap.clear();
    
    //read info from stable file
    string testfile = testDir+"test.names";
    util.readNames(testfile, nameMap);
}
/***********************************************************************/
void TestDataSet::fillGroup() {
    if (gMap != NULL) { delete gMap; gMap = NULL; }
    
    //read info from stable file
    string testfile = testDir+"test.groups";
    
    gMap = new GroupMap();
    gMap->readMap(testfile);
}
/***********************************************************************/
void TestDataSet::fillLookup() {
    for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  lookup[i] = NULL; }
    lookup.clear();
    
    //read info from stable file
    string testfile = testDir+"test.opti_mcc.shared";

    InputData input(testfile, "sharedfile", nullVector);
    SharedRAbundVectors* shared = input.getSharedRAbundVectors();
    lookup = shared->getSharedRAbundVectors();
}
/***********************************************************************/


