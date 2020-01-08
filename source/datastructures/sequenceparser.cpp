/*
 *  sequenceParser.cpp
 *  Mothur
 *
 *  Created by westcott on 9/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "sequenceparser.h"

/************************************************************/
SequenceParser::SequenceParser(string groupFile, string fastaFile, string nameFile, vector<string> groupsSelected) {
	try {
		
		m = MothurOut::getInstance();
        hasName = true;
		
		//read group file
        GroupMap groupMap;
		int error = groupMap.readMap(groupFile, groupsSelected); //only store info for groups selected
		
		if (error == 1) { m->setControl_pressed(true); }
		
		//initialize maps
        namesOfGroups = groupMap.getNamesOfGroups();
        
        //run splitGroups command to parse files
        string inputString = "";
        if (groupsSelected.size() != 0) { sort(groupsSelected.begin(), groupsSelected.end());   }
        else                            { groupsSelected = namesOfGroups;                       }
        
        inputString += "processors=1, groups=" + util.getStringFromVector(groupsSelected, "-"); //split.groups is paraplellized, we don't want the thread spinning up threads.
        inputString += ", fasta=" + fastaFile;
        inputString += ", name=" + nameFile;
        inputString += ", group=" + groupFile;
        
        m->mothurOut("\n/******************************************/\n");
        m->mothurOut("Running command: split.groups(" + inputString + ")\n");
        
        Command* splitCommand = new SplitGroupCommand(inputString);
        splitCommand->execute();
        
        //type -> files in groups order. fasta -> vector<string>. fastaFileForGroup1 stored in filenames["fasta"][1]
        map<string, vector<string> > filenames = splitCommand->getOutputFiles();
        
        delete splitCommand;
        m->mothurOut("/******************************************/\n");
        
        vector<string> parsedFastaFiles = filenames["fasta"]; //sorted in groups order
        vector<string> parsedNameFiles = filenames["name"]; //sorted in groups order
        vector<string> parsedGroupFiles = filenames["group"]; //sorted in groups order
        
        if (parsedNameFiles.size() != groupsSelected.size()) { cout << "should never get here, quitting\n\n"; m->setControl_pressed(true);  }
        
        for (int i = 0; i < groupsSelected.size(); i++) {
            vector<string> thisSamplesFiles;
            thisSamplesFiles.push_back(parsedFastaFiles[i]);
            thisSamplesFiles.push_back(parsedNameFiles[i]);
            thisSamplesFiles.push_back(parsedGroupFiles[i]);
            groupToFiles[groupsSelected[i]] = thisSamplesFiles;
        }
        
        //reset current files changed by split.groups
        CurrentFile* current; current = CurrentFile::getInstance();
        current->setNameFile(nameFile);
        current->setFastaFile(fastaFile);
        current->setGroupFile(groupFile);
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "SequenceParser");
		exit(1);
	}
}
/************************************************************/
//leaves all seqs map blank to be filled when asked for
SequenceParser::SequenceParser(string groupFile, string fastaFile, vector<string> groupsSelected) {
	try {
		
        m = MothurOut::getInstance();
        hasName = false;
        
        //read group file
        GroupMap groupMap;
        int error = groupMap.readMap(groupFile, groupsSelected); //only store info for groups selected
        
        if (error == 1) { m->setControl_pressed(true); }
        
        //initialize maps
        namesOfGroups = groupMap.getNamesOfGroups();
        
        //run splitGroups command to parse files
        string inputString = "";
        if (groupsSelected.size() != 0) {
            sort(groupsSelected.begin(), groupsSelected.end());
            inputString += "groups=" + util.getStringFromVector(groupsSelected, "-");
        }else {
            groupsSelected = namesOfGroups;
        }
        
        inputString += ", fasta=" + fastaFile;
        inputString += ", group=" + groupFile;
        
        m->mothurOut("\n/******************************************/\n");
        m->mothurOut("Running command: split.groups(" + inputString + ")\n");
        
        Command* splitCommand = new SplitGroupCommand(inputString);
        splitCommand->execute();
        
        //type -> files in groups order. fasta -> vector<string>. fastaFileForGroup1 stored in filenames["fasta"][1]
        map<string, vector<string> > filenames = splitCommand->getOutputFiles();
        
        delete splitCommand;
        m->mothurOut("/******************************************/\n");
        
        vector<string> parsedFastaFiles = filenames["fasta"]; //sorted in groups order
        vector<string> parsedGroupFiles = filenames["group"]; //sorted in groups order
        
        if (parsedFastaFiles.size() != groupsSelected.size()) { cout << "should never get here, quitting\n\n"; m->setControl_pressed(true);  }
        
        for (int i = 0; i < groupsSelected.size(); i++) {
            vector<string> thisSamplesFiles;
            thisSamplesFiles.push_back(parsedFastaFiles[i]);
            thisSamplesFiles.push_back(parsedGroupFiles[i]);
            groupToFiles[groupsSelected[i]] = thisSamplesFiles;
        }
        //reset current files changed by split.groups
        CurrentFile* current; current = CurrentFile::getInstance();
        current->setFastaFile(fastaFile);
        current->setGroupFile(groupFile);
		
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "SequenceParser");
		exit(1);
	}
}
/************************************************************/
SequenceParser::~SequenceParser(){ }
/************************************************************/
int SequenceParser::getNumGroups(){ return namesOfGroups.size(); }
/************************************************************/
vector<string> SequenceParser::getNamesOfGroups(){ return namesOfGroups; }
/************************************************************/
vector<string> SequenceParser::getFiles(string group){
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
        m->errorOut(e, "SequenceParser", "getFiles");
        exit(1);
    }
}
/************************************************************/


