//
//  sequencecountparser.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/7/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sequencecountparser.h"
#include "sharedutilities.h"

/************************************************************/
SequenceCountParser::SequenceCountParser(string countfile, string fastafile, vector<string> groupsSelected) {
	try {
		
		m = MothurOut::getInstance();
		
		//read count file
		CountTable countTable;
		countTable.readTable(countfile, true, false);
        vector<string> allNames = countTable.getNamesOfGroups();
		
		//initialize maps
        if (groupsSelected.size() == 0) { //select all
            namesOfGroups = countTable.getNamesOfGroups();
            for (int i = 0; i < allNames.size(); i++) { indexes.push_back(i); }
        }else{
            SharedUtil util;  util.setGroups(groupsSelected, namesOfGroups);
            namesOfGroups = groupsSelected;
            map<string, int> temp;
            for (int i = 0; i < allNames.size(); i++) {
                for (int j = 0; j < groupsSelected.size(); j++) {
                    if (allNames[i] == groupsSelected[j]) {
                        temp[groupsSelected[j]] = i;
                    }
                }
            }
            for (map<string, int>::iterator it = temp.begin(); it != temp.end(); it++) {
                indexes.push_back(it->second);
            }
        }
        
		for (int i = 0; i < namesOfGroups.size(); i++) {
			vector<Sequence> temp;
			map<string, int> tempMap;
			seqs[namesOfGroups[i]] = temp;
			countTablePerGroup[namesOfGroups[i]] = tempMap;
		}
        
        readFasta(fastafile, countTable);
    }
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "SequenceCountParser");
		exit(1);
	}
}
/************************************************************/
SequenceCountParser::SequenceCountParser(string fastafile, CountTable& countTable, vector<string> groupsSelected) {
	try {
		
		m = MothurOut::getInstance();
				
		//initialize maps
        if (countTable.hasGroupInfo()) {
            vector<string> allNames = countTable.getNamesOfGroups();
            
            //initialize maps
            if (groupsSelected.size() == 0) { //select all
                namesOfGroups = countTable.getNamesOfGroups();
                for (int i = 0; i < allNames.size(); i++) { indexes.push_back(i); }
            }else{
                SharedUtil util;  util.setGroups(groupsSelected, namesOfGroups);
                namesOfGroups = groupsSelected;
                map<string, int> temp;
                for (int i = 0; i < allNames.size(); i++) {
                    for (int j = 0; j < groupsSelected.size(); j++) {
                        if (allNames[i] == groupsSelected[j]) {
                            temp[groupsSelected[j]] = i;
                        }
                    }
                }
                for (map<string, int>::iterator it = temp.begin(); it != temp.end(); it++) {
                    indexes.push_back(it->second);
                }
            }

            for (int i = 0; i < namesOfGroups.size(); i++) {
                vector<Sequence> temp;
                map<string, int> tempMap;
                seqs[namesOfGroups[i]] = temp;
                countTablePerGroup[namesOfGroups[i]] = tempMap;
            }
            
            readFasta(fastafile, countTable);
            
        }else {  m->control_pressed = true;  m->mothurOut("[ERROR]: cannot parse fasta file by group with a count table that does not include group data, please correct.\n"); }
        
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "SequenceCountParser");
		exit(1);
	}
}
/************************************************************/
int SequenceCountParser::readFasta(string fastafile, CountTable& countTable) {
    try {
        
        
        ifstream in;
        m->openInputFile(fastafile, in);
        
        int fastaCount = 0;
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            Sequence seq(in); m->gobble(in);
            fastaCount++;
            if (m->debug) { if((fastaCount) % 1000 == 0){	m->mothurOut("[DEBUG]: reading seq " + toString(fastaCount) + "\n.");	} }
            
            if (seq.getName() != "") {
                
                //allSeqsMap[seq.getName()] = seq.getName();
                vector<int> groupCounts = countTable.getGroupCounts(seq.getName());
                
                for (int i = 0; i < namesOfGroups.size(); i++) {
                    if (groupCounts[indexes[i]] != 0) {
                        seqs[namesOfGroups[i]].push_back(seq);
                        countTablePerGroup[namesOfGroups[i]][seq.getName()] = groupCounts[i];
                    }
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
		map<string, vector<Sequence> >::iterator it;
		int num = 0;
		
		it = seqs.find(g);
		if(it == seqs.end()) {
			m->mothurOut("[ERROR]: " + g + " is not a valid group, please correct."); m->mothurOutEndLine();
		}else {
			num = (it->second).size();
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
        map<string, vector<Sequence> >::iterator it;
        
        for (it = seqs.begin(); it != seqs.end(); it++) {
            vector<Sequence> temp = it->second;
            for (int i = 0; i < temp.size(); i++) {  allSeqsMap[temp[i].getName()] = temp[i].getName(); }
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
		map<string, vector<Sequence> >::iterator it;
		vector<Sequence> seqForThisGroup;
		
		it = seqs.find(g);
		if(it == seqs.end()) {
			m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			seqForThisGroup = it->second;
            if (m->debug) {  m->mothurOut("[DEBUG]: group " + g + " fasta file has " + toString(seqForThisGroup.size()) + " sequences.");  }
		}
		
		return seqForThisGroup; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceCountParser", "getSeqs");
		exit(1);
	}
}
/************************************************************/
int SequenceCountParser::getSeqs(string g, string filename, string tag, string tag2, bool uchimeFormat=false){
	try {
		map<string, vector<Sequence> >::iterator it;
		vector<Sequence> seqForThisGroup;
		vector<seqPriorityNode> nameVector;
		
		it = seqs.find(g);
		if(it == seqs.end()) {
			m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			
			ofstream out;
			m->openOutputFile(filename, out);
			
			seqForThisGroup = it->second;
			
			if (uchimeFormat) {
				// format should look like 
				//>seqName /ab=numRedundantSeqs/
				//sequence
				
				map<string, int> countForThisGroup = getCountTable(g);
				map<string, int>::iterator itCount;
				int error = 0;
				
				for (int i = 0; i < seqForThisGroup.size(); i++) {
					itCount = countForThisGroup.find(seqForThisGroup[i].getName());
					
					if (itCount == countForThisGroup.end()){
						error = 1;
						m->mothurOut("[ERROR]: " + seqForThisGroup[i].getName() + " is in your fastafile, but is not in your count file, please correct."); m->mothurOutEndLine();
					}else {
                        seqPriorityNode temp(itCount->second, seqForThisGroup[i].getUnaligned(), seqForThisGroup[i].getName());
						nameVector.push_back(temp);
					}
				}
				
				if (error == 1) { out.close(); m->mothurRemove(filename); return 1; }
				
				//sort by num represented
				sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);
                
				//print new file in order of
				for (int i = 0; i < nameVector.size(); i++) {
					
					if(m->control_pressed) { out.close(); m->mothurRemove(filename); return 1; }
					
					out << ">" << nameVector[i].name << tag << nameVector[i].numIdentical << tag2 << endl << nameVector[i].seq << endl; //
				}
				
			}else { 
                //m->mothurOut("Group " + g +  " contains " + toString(seqForThisGroup.size()) + " unique seqs.\n");
				for (int i = 0; i < seqForThisGroup.size(); i++) {
					
					if(m->control_pressed) { out.close(); m->mothurRemove(filename); return 1; }
					
					seqForThisGroup[i].printSequence(out);	
				}
			}
			out.close();
		}
		
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
		map<string, map<string, int> >::iterator it;
		map<string, int> countForThisGroup;
		
		it = countTablePerGroup.find(g);
		if(it == countTablePerGroup.end()) {
			m->mothurOut("[ERROR]: No countTable available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			countForThisGroup = it->second;
            if (m->debug) {  m->mothurOut("[DEBUG]: group " + g + " count file has " + toString(countForThisGroup.size()) + " unique sequences.");  }
		}
		
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
		map<string, map<string, int> >::iterator it;
		map<string, int> countForThisGroup;
		
		it = countTablePerGroup.find(g);
		if(it == countTablePerGroup.end()) {
			m->mothurOut("[ERROR]: No countTable available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			countForThisGroup = it->second;
			
			ofstream out;
			m->openOutputFile(filename, out);
            out << "Representative_Sequence\ttotal\t" << g << endl;
            
			for (map<string, int>::iterator itFile = countForThisGroup.begin(); itFile != countForThisGroup.end(); itFile++) {
				
				if(m->control_pressed) { out.close(); m->mothurRemove(filename); return 1; }
				
				out << itFile->first << '\t' << itFile->second << '\t' << itFile->second << endl;
			}
			
			out.close();
		}
		
		return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getCountTable");
		exit(1);
	}
}
/************************************************************/



