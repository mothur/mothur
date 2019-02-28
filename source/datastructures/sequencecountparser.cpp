//
//  sequencecountparser.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/7/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sequencecountparser.h"


/************************************************************/
SequenceCountParser::SequenceCountParser(string countfile, string fastafile, vector<string> groupsSelected) {
	try {
        
		m = MothurOut::getInstance();
		
		//read count file
		countTable.readTable(countfile, true, false, groupsSelected);
        namesOfGroups = countTable.getNamesOfGroups();

        groupToSeqs.resize(namesOfGroups.size()); //allocate space for samples
        
        for (int i = 0; i < namesOfGroups.size(); i++) { groupIndexMap[namesOfGroups[i]] = i; }
        
        readFasta(fastafile);
    }
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "SequenceCountParser");
		exit(1);
	}
}
/************************************************************/
SequenceCountParser::SequenceCountParser(string fastafile, CountTable ct, vector<string> groupsSelected) {
    try {
        
        m = MothurOut::getInstance();
        
        //read count file
        countTable.copy(&ct);
        vector<string> allNames = countTable.getNamesOfGroups();
        
        //initialize maps
        if (groupsSelected.size() == 0) { //select all
            namesOfGroups = countTable.getNamesOfGroups();
        }else{
            namesOfGroups = groupsSelected;
            //remove groups from count table that we don't want
            for (int i = 0; i < allNames.size(); i++) {
                if (!util.inUsersGroups(allNames[i], groupsSelected)) { //we want to remove unwanted groups from memory
                    countTable.removeGroup(allNames[i]);
                }
            }
        }
        
        groupToSeqs.resize(namesOfGroups.size()); //allocate space for samples
        
        for (int i = 0; i < namesOfGroups.size(); i++) { groupIndexMap[namesOfGroups[i]] = i; }
        
        readFasta(fastafile);
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceCountParser", "SequenceCountParser");
        exit(1);
    }
}
/************************************************************/
int SequenceCountParser::readFasta(string fastafile) {
    try {
        ifstream in;
        Utils util; util.openInputFile(fastafile, in);
        
        int count = 0;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); util.gobble(in);
            
            if (seq.getName() != "") {
                
                if (countTable.inTable(seq.getName())) { //only save sequences from groups we want
                    vector<int> groupCounts = countTable.getGroupCounts(seq.getName());
                    
                    seqs.push_back(seq);
                    
                    for (int i = 0; i < namesOfGroups.size(); i++) {
                        if (groupCounts[i] != 0) { groupToSeqs[i].push_back(count); }
                    }
                    
                    count++;
                }
            }
        }
        in.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceCountParser", "readFasta");
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
int SequenceCountParser::getNumSeqs(string g){ 
	try {
		map<string, int>::iterator it;
		int num = 0;
		
		it = groupIndexMap.find(g);
		if(it == groupIndexMap.end()) {
			m->mothurOut("[ERROR]: " + g + " is not a valid group, please correct.\n"); 
		}else {
			num = groupToSeqs[it->second].size();
		}
		
		return num; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "getNumSeqs");
		exit(1);
	}
}
/************************************************************/
map<string, string> SequenceCountParser::getAllSeqsMap(){
    try {
        map<string, string > allSeqsMap;
        
        for (int i = 0; i < seqs.size(); i++) {
            allSeqsMap[seqs[i].getName()] = seqs[i].getName();
        }
        
        return allSeqsMap;
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceCountParser", "getAllSeqsMap");
        exit(1);
    }
}
/************************************************************/
vector<Sequence> SequenceCountParser::getSeqs(string g){ 
	try {
		
		vector<Sequence> seqForThisGroup;
        map<string, int>::iterator it;
        
        it = groupIndexMap.find(g);
        if(it == groupIndexMap.end()) {
            m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct.\n");
        }else {
            for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                seqForThisGroup.push_back(seqs[groupToSeqs[it->second][i]]);
            }
        }
		
		return seqForThisGroup; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "getSeqs");
		exit(1);
	}
}
/************************************************************/
bool SequenceCountParser::fillWeighted(vector< seqPNode* >& seqForThisGroup, string group, int& length){
    try {
        
        set<int> lengths;
        
        map<string, int>::iterator it;
        
        it = groupIndexMap.find(group);
        if(it == groupIndexMap.end()) {
            m->mothurOut("[ERROR]: No sequences available for group " + group + ", please correct.\n"); return false;
        }else {
            for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                
                Sequence thisSeq = seqs[groupToSeqs[it->second][i]];
                int numReps = countTable.getGroupCount(thisSeq.getName(), group);
                
                set<int> clusteredIndexes; //don't use indexes for precluster with count file
                seqPNode* tempNode = new seqPNode(thisSeq.getName(), thisSeq.getAligned(), numReps, clusteredIndexes);
                seqForThisGroup.push_back(tempNode);
                
                lengths.insert(thisSeq.getAligned().length());
            }
        }

        length = *(lengths.begin());
        
        if (lengths.size() > 1) { return false; } //unaligned
        else if (lengths.size() == 1) {  return true; } //aligned
        
        return true;
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceCountParser", "fillWeighted");
        exit(1);
    }
}
/************************************************************/
int SequenceCountParser::getSeqs(string g, string filename, string tag, string tag2, long long& numSeqs, bool uchimeFormat=false){
	try {
        vector<seqPriorityNode> nameVector;
        
        ofstream out;
        Utils util; util.openOutputFile(filename, out);
        
        vector<Sequence> seqForThisGroup = getSeqs(g);
        
        numSeqs = seqForThisGroup.size();
        
        if (uchimeFormat) {
            // format should look like
            //>seqName /ab=numRedundantSeqs/
            //sequence
            
            map<string, int> countForThisGroup = getCountTable(g);
            map<string, int>::iterator itCount;
            bool error = false;
            
            for (int i = 0; i < seqForThisGroup.size(); i++) {
                itCount = countForThisGroup.find(seqForThisGroup[i].getName());
                
                if (itCount == countForThisGroup.end()){
                    error = true;
                    m->mothurOut("[ERROR]: " + seqForThisGroup[i].getName() + " is in your fastafile, but is not in your count file, please correct."); m->mothurOutEndLine();
                }else {
                    seqPriorityNode temp(itCount->second, seqForThisGroup[i].getUnaligned(), seqForThisGroup[i].getName());
                    nameVector.push_back(temp);
                }
            }
            
            if (error) { out.close(); util.mothurRemove(filename); return 1; }
            
            //sort by num represented
            sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);
            
            //print new file in order of
            for (int i = 0; i < nameVector.size(); i++) {
                
                if(m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 1; }
                
                out << ">" << nameVector[i].name << tag << nameVector[i].numIdentical << tag2 << endl << nameVector[i].seq << endl; //
            }
            
        }else {
            //m->mothurOut("Group " + g +  " contains " + toString(seqForThisGroup.size()) + " unique seqs.\n");
            for (int i = 0; i < seqForThisGroup.size(); i++) {
                
                if(m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 1; }
                
                seqForThisGroup[i].printSequence(out);	
            }
        }
        out.close();
        
        
		return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "getSeqs");
		exit(1);
	}
}

/************************************************************/
map<string, int> SequenceCountParser::getCountTable(string g){ 
	try {
		map<string, int> countForThisGroup = countTable.getNameMap(g);
		
		return countForThisGroup; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "getCountTable");
		exit(1);
	}
}
/************************************************************/
int SequenceCountParser::getCountTable(string g, string filename){ 
	try {
        Utils util;
        map<string, int> countForThisGroup = getCountTable(g);
        
        ofstream out;
        util.openOutputFile(filename, out);
        out << "Representative_Sequence\ttotal\t" << g << endl;
        
        for (map<string, int>::iterator itFile = countForThisGroup.begin(); itFile != countForThisGroup.end(); itFile++) {
            
            if(m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 1; }
            
            out << itFile->first << '\t' << itFile->second << '\t' << itFile->second << endl;
        }
        
        out.close();
        
		return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getCountTable");
		exit(1);
	}
}
/************************************************************/



