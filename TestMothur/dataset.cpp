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
        m->splitAtChar(secondCol, names, ',');
        
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

vector<string> TestDataSet::getSubsetFNGFiles(int numSeqs) {
    fillSeqs();
    vector<Sequence> subsetSeqs;
    for (int i = 0; i < numSeqs; i++) { subsetSeqs.push_back(seqs[i]); }
    seqs.clear();
    
    fillNames();
    fillGroup();
    ofstream out, out2, out3;
    m->openOutputFile("tempSeqs.txt", out); m->openOutputFile("tempNames.txt", out2); m->openOutputFile("tempGroup.txt", out3);
    for (int i = 0; i < subsetSeqs.size(); i++) {
        subsetSeqs[i].printSequence(out);
        out2 << subsetSeqs[i].getName() << '\t' << nameMap[subsetSeqs[i].getName()] << '\n';
        out3 << subsetSeqs[i].getName() << '\t' << gMap->getGroup(subsetSeqs[i].getName()) << '\n';
    }
    nameMap.clear();
    delete gMap; gMap = NULL;
    
    vector<string> filenames; filenames.push_back("tempSeqs.txt"); filenames.push_back("tempNames.txt"); filenames.push_back("tempGroup.txt");
    
    return filenames;
}
/***********************************************************************/
void TestDataSet::fillSeqs() {
    seqs.clear();
    
    //read info from stable file
    //string testfile = m->getTestFilePath() + "testFile.fasta";
     string testfile ="/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/testFile.fasta";
    
    ifstream in;
    m->openInputFile(testfile, in);
    
    while (!in.eof()) {
        if (m->control_pressed) { break; }
        
        Sequence read(in); m->gobble(in);
        seqs.push_back(read);
    }
    in.close();
}
/***********************************************************************/
void TestDataSet::fillNames() {
    nameMap.clear();
    
    //read info from stable file
    string testfile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/testFile.names";
    m->readNames(testfile, nameMap);
}
/***********************************************************************/
void TestDataSet::fillGroup() {
    if (gMap != NULL) { delete gMap; gMap = NULL; }
    
    //read info from stable file
    string testfile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/testFile.groups";
    
    gMap = new GroupMap();
    gMap->readMap(testfile);
}
/***********************************************************************/
void TestDataSet::fillLookup() {
    for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  lookup[i] = NULL; }
    lookup.clear();
    
    //read info from stable file
    string testfile = "/Users/sarahwestcott/Desktop/mothur/TestMothur/TestFiles/testFile.opti_mcc.shared";

    InputData input(testfile, "sharedfile");
    lookup = input.getSharedRAbundVectors();
}
/***********************************************************************/


