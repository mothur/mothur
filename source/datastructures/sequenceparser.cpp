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
        
        inputString += "groups=" + util.getStringFromVector(groupsSelected, "-");
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
/************************************************************
int SequenceParser::getSeqs(string g, string filename, string tag, string tag2, long long& numSeqs, string optionalGroupFileName, bool uchimeFormat=false){
	try {
		map<string, int>::iterator it;
		vector<Sequence> seqForThisGroup;
		vector<seqPriorityNode> nameVector;
		
        it = groupIndexMap.find(g);
        if(it == groupIndexMap.end()) {
            m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct.\n");
        }else {

			ofstream out;
			util.openOutputFile(filename, out);
			
			seqForThisGroup = getSeqs(g);
            
            numSeqs = seqForThisGroup.size();

			if (uchimeFormat) {
				// format should look like 
				//>seqName /ab=numRedundantSeqs/
				//sequence
				
				map<string, string> nameMapForThisGroup = getNameMap(g);
				map<string, string>::iterator itNameMap;
				int error = 0;
				
				for (int i = 0; i < seqForThisGroup.size(); i++) {
					itNameMap = nameMapForThisGroup.find(seqForThisGroup[i].getName());
					
					if (itNameMap == nameMapForThisGroup.end()){
						error = 1;
						m->mothurOut("[ERROR]: " + seqForThisGroup[i].getName() + " is in your fastafile, but is not in your namesfile, please correct."); m->mothurOutEndLine();
					}else {
						int num = util.getNumNames(itNameMap->second);
						
						seqPriorityNode temp(num, seqForThisGroup[i].getUnaligned(), seqForThisGroup[i].getName());
						nameVector.push_back(temp);
					}
				}
				
				if (error == 1) { out.close(); util.mothurRemove(filename); return 1; }
				
				//sort by num represented
				sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);

				//print new file in order of
				for (int i = 0; i < nameVector.size(); i++) {
					
					if(m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 1; }
					
					out << ">" <<  nameVector[i].name << tag << nameVector[i].numIdentical << tag2 << endl << nameVector[i].seq << endl; //
				}
				
			}else {
                
                ofstream outGroup;
                if (optionalGroupFileName != "") { util.openOutputFile(optionalGroupFileName, outGroup); }
                
				for (int i = 0; i < seqForThisGroup.size(); i++) {
					
					if(m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 1; }
					
					seqForThisGroup[i].printSequence(out);
                    
                    if (optionalGroupFileName != "") {  outGroup << seqForThisGroup[i].getName() << '\t' << g << endl;  }
				}
                
                if (optionalGroupFileName != "") {  outGroup.close();  }
			}
			out.close();
            
		}
		
		return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getSeqs");
		exit(1);
	}
}
/************************************************************/



