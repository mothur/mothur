//
//  sequencecountparser.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/7/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sequencecountparser.h"
#include "splitgroupscommand.h"

/************************************************************/
SequenceCountParser::SequenceCountParser(string countfile, string fastafile, vector<string> groupsSelected) {
	try {
        
		m = MothurOut::getInstance();
    
        //run splitGroups command to parse files
        string inputString = "";
        if (groupsSelected.size() == 0) {
            CountTable ct; ct.testGroups(countfile, groupsSelected); //fills groupsSelected with groups in count table
        }
        
        m->mothurOut("\n/******************************************/\n");
        m->mothurOut("Splitting by sample: \n");
        
        time_t start = time(nullptr);
        
        SplitGroupCommand* splitCommand = new SplitGroupCommand(groupsSelected, fastafile, countfile, "");

        //type -> files in groups order. fasta -> vector<string>. fastaFileForGroup1 stored in filenames["fasta"][1]
        map<string, vector<string> > filenames = splitCommand->getOutputFiles();
        
        delete splitCommand;
        
        m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " seconds to split the dataset by sample.\n");
        
        m->mothurOut("/******************************************/\n");
        
        vector<string> parsedFastaFiles = filenames["fasta"]; //sorted in groups order
        vector<string> parsedCountFiles = filenames["count"]; //sorted in groups order
        
        if (parsedCountFiles.size() != groupsSelected.size()) { cout << "should never get here, quitting\n\n"; m->setControl_pressed(true);  }
        
        namesOfGroups = groupsSelected;
        
        for (int i = 0; i < groupsSelected.size(); i++) {
            vector<string> thisSamplesFiles;
            thisSamplesFiles.push_back(parsedFastaFiles[i]);
            thisSamplesFiles.push_back(parsedCountFiles[i]);
            groupToFiles[groupsSelected[i]] = thisSamplesFiles;
        }
        
        //reset current files changed by split.groups
        CurrentFile* current; current = CurrentFile::getInstance();
        current->setCountFile(countfile);
        current->setFastaFile(fastafile);
        
    }
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "SequenceCountParser");
		exit(1);
	}
}
/************************************************************/
SequenceCountParser::~SequenceCountParser(){  }
/************************************************************/
int SequenceCountParser::getNumGroups(){ return namesOfGroups.size(); }
/************************************************************/
vector<string> SequenceCountParser::getNamesOfGroups(){ return namesOfGroups; }
/************************************************************/
vector<string> SequenceCountParser::getFiles(string group){
    try {
        map<string, vector<string> >::iterator it;
        
        it = groupToFiles.find(group);
        if (it != groupToFiles.end()) {
            return it->second;
        }else {
            m->mothurOut("[ERROR]: cannot find files for group " + group + ", quitting.\n"); m->setControl_pressed(true);
        }
        
        return nullVector;
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceCountParser", "getFiles");
        exit(1);
    }
}
/************************************************************/


