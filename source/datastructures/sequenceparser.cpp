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
		int error = groupMap.readMap(groupFile, groupsSelected); //only store info for groups selected
		
		if (error == 1) { m->setControl_pressed(true); }
		
		//initialize maps
        namesOfGroups = groupMap.getNamesOfGroups();
        
        groupToSeqs.resize(namesOfGroups.size()); //allocate space for samples
        
        for (int i = 0; i < namesOfGroups.size(); i++) { groupIndexMap[namesOfGroups[i]] = i; }
        
        map<string, string>::iterator it;
        vector<string> input = groupMap.getNamesSeqs();
        set<string> namesToInclude = util.mothurConvert(input); input.clear();
        util.readNames(nameFile, nameMap, namesToInclude); //only reads names included in group map
        input.clear();
        
        int count = 0;
        
		//read fasta file making sure each sequence is in the group file
		ifstream in;
		util.openInputFile(fastaFile, in);
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); util.gobble(in);
            
            if (seq.getName() != "") {
                
                it = nameMap.find(seq.getName());
                
                if (it != nameMap.end()) { //in namefile and users selected groups
                    
                    seqs.push_back(seq);
                    
                    vector<string> names;
                    string secondCol = it->second;
                    util.splitAtChar(secondCol, names, ',');
                    
                    //fills allSeqsMap
                    for (int i = 0; i < names.size(); i++) { allSeqsMap[names[i]] = names[0]; }
                    
                    //fill groupsToSeqs
                    vector<string> representedGroups = groupMap.getGroups(names);
                    for (int i = 0; i < representedGroups.size(); i++) {
                        map<string, int>::iterator itGroupIndex = groupIndexMap.find(representedGroups[i]); //find index of group in groupsToSeqs
                        
                        if (itGroupIndex != groupIndexMap.end()) {
                            groupToSeqs[itGroupIndex->second].push_back(count);
                        }
                    }
                    
                    count++;
                }//else ignore seq, its not from groups we want
            }
        }
		in.close();
				 
		if (error == 1) { m->setControl_pressed(true); }
		
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
        int error = groupMap.readMap(groupFile, groupsSelected); //only store info for groups selected
        
        if (error == 1) { m->setControl_pressed(true); }
        
        //initialize maps
        namesOfGroups = groupMap.getNamesOfGroups();
        
        groupToSeqs.resize(namesOfGroups.size()); //allocate space for samples
        
        for (int i = 0; i < namesOfGroups.size(); i++) { groupIndexMap[namesOfGroups[i]] = i; }

		//read fasta file making sure each sequence is in the group file
		ifstream in;
		util.openInputFile(fastaFile, in);
		
        int count = 0;
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { break; }
			
			Sequence seq(in); util.gobble(in);
			
			if (seq.getName() != "") {
                
                seqs.push_back(seq);
                
                //fill groupsToSeqs
                string group = groupMap.getGroup(seq.getName());
                
                map<string, int>::iterator itGroupIndex = groupIndexMap.find(group); //find index of group in groupsToSeqs
                
                if (itGroupIndex != groupIndexMap.end()) {
                    groupToSeqs[itGroupIndex->second].push_back(count);
                }
                
                count++;

			}
		}
		in.close();
		
		if (error == 1) { m->setControl_pressed(true); }
		
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
int SequenceParser::getNumSeqs(string g){ 
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
		m->errorOut(e, "SequenceParser", "getNumSeqs");
		exit(1);
	}
}
/************************************************************/
vector<Sequence> SequenceParser::getSeqs(string g){ 
	try {
        vector<Sequence> seqForThisGroup;
        map<string, int>::iterator it;
        
        it = groupIndexMap.find(g);
        if(it == groupIndexMap.end()) {
            m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct.\n");
        }else {
            if (hasName) {
                for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                    string uniqueSeqName = seqs[groupToSeqs[it->second][i]].getName();
                    string aligned = seqs[groupToSeqs[it->second][i]].getAligned();
                    
                    string dupNameForThisGroup = ""; //we want to set this to the first seq name we find in the dups names from group g
                    
                    map<string, string>::iterator it = nameMap.find(uniqueSeqName);
                    if (it != nameMap.end()) {
                        string dups = it->second;
                        vector<string> dupsNames; util.splitAtComma(dups, dupsNames);
                        
                        for (int j = 0; j < dupsNames.size(); j++) {
                            string dupsGroup = groupMap.getGroup(dupsNames[j]);
                            if (dupsGroup == g) {
                                dupNameForThisGroup = dupsNames[j];
                                break;
                            }
                        }
                    }
                    
                    if (dupNameForThisGroup != "") {
                        Sequence thisDupSeq(dupNameForThisGroup, aligned);
                        seqForThisGroup.push_back(thisDupSeq);
                    }else { m->mothurOut("[ERROR]: should never get here\n");  m->setControl_pressed(true); }
                }
            }else { //seq names are unique
                for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                    seqForThisGroup.push_back(seqs[groupToSeqs[it->second][i]]);
                }
            }
        }
        
        return seqForThisGroup;
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getSeqs");
		exit(1);
	}
}
/************************************************************/
bool SequenceParser::fillWeighted(vector< seqPNode* >& seqForThisGroup, string group, int& length){
    try {
        set<int> lengths;
        
        map<string, int>::iterator it;
        
        it = groupIndexMap.find(group);
        if(it == groupIndexMap.end()) {
            m->mothurOut("[ERROR]: No sequences available for group " + group + ", please correct.\n"); return false;
        }else {
            for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                
                Sequence thisSeq = seqs[groupToSeqs[it->second][i]];
                int numReps = 1;
                
                string uniqueSeqNameForGroup = thisSeq.getName(); //if the representative sequence for this read is not from the group we are looking for we need to change the name
                if (hasName) {
                    //find your nameFile dups
                    map<string, string>::iterator it = nameMap.find(thisSeq.getName());
                    if (it != nameMap.end()) {
                        string dups = it->second;
                        numReps = groupMap.getNumSeqs(dups, group);
                        
                        string d = it->second;
                        vector<string> dupsNames; util.splitAtComma(d, dupsNames);
                        for ( int j = 0; j < dupsNames.size(); j++) {
                            if (groupMap.getGroup(dupsNames[j]) == group) {
                                uniqueSeqNameForGroup = dupsNames[j]; break;
                            }
                        }
                    }
                }
                
                set<int> clusteredIndexes; clusteredIndexes.insert(seqForThisGroup.size());
                seqPNode* tempNode = new seqPNode(uniqueSeqNameForGroup, thisSeq.getAligned(), numReps, clusteredIndexes);
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
        m->errorOut(e, "SequenceParser", "fillWeighted");
        exit(1);
    }
}
/************************************************************/
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
map<string, string> SequenceParser::getAllSeqsMap(){
    try {
        map<string, string> uniqueAllSeqs;
        
        if (hasName) { return allSeqsMap; }
        else { //create unique allSeqsMap
            for (int i = 0; i < seqs.size(); i++) {
                uniqueAllSeqs[seqs[i].getName()] = seqs[i].getName();
            }
        }
        return uniqueAllSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceParser", "getAllSeqsMap");
        exit(1);
    }
}
/************************************************************/
map<string, string> SequenceParser::getNameMap(string g){ 
	try {
		map<string, string> nameMapForThisGroup;
        map<string, int>::iterator it;
        
        it = groupIndexMap.find(g);
        if(it == groupIndexMap.end()) {
            m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct.\n"); return nameMapForThisGroup;
        }else {

            if (hasName) {
                for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                    string uniqueSeqName = seqs[groupToSeqs[it->second][i]].getName();
                    string secondCol = "";
                    
                    map<string, string>::iterator it = nameMap.find(uniqueSeqName);
                    if (it != nameMap.end()) {
                        string dups = it->second; //contains all groups
                        vector<string> dupsNames; util.splitAtComma(dups, dupsNames);
                        
                        for (int j = 0; j < dupsNames.size(); j++) {
                            string dupsGroup = groupMap.getGroup(dupsNames[j]);
                            if (dupsGroup == g) {
                                secondCol += dupsNames[j] +",";
                            }
                        }
                    }
                    
                    if (secondCol != "") {
                        //remove last comma
                        secondCol = secondCol.substr(0,secondCol.length()-1);
                        int pos = secondCol.find_first_of(',');
                        string firstCol = secondCol;
                        if (pos != string::npos) {
                            firstCol = secondCol.substr(0, pos);
                        }
                        nameMapForThisGroup[firstCol] = secondCol;
                        
                    }else { m->mothurOut("[ERROR]: should never get here\n"); m->setControl_pressed(true); }
                }
            }else { //seq names are unique
                for (int i = 0; i < groupToSeqs[it->second].size(); i++) {
                    nameMapForThisGroup[seqs[groupToSeqs[it->second][i]].getName()] = seqs[groupToSeqs[it->second][i]].getName();
                }
            }

		}
		
		return nameMapForThisGroup; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getNameMap");
		exit(1);
	}
}
/************************************************************/
int SequenceParser::getNameMap(string g, string filename){ 
	try {
		
        map<string, string> nameMapForThisGroup = getNameMap(g);

        ofstream out;
        util.openOutputFile(filename, out);
        
        for (map<string, string>::iterator itFile = nameMapForThisGroup.begin(); itFile != nameMapForThisGroup.end(); itFile++) {
            
            if(m->getControl_pressed()) { out.close(); util.mothurRemove(filename); return 1; }
            
            out << itFile->first << '\t' << itFile->second << endl;
        }
        
        out.close();
		
		return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getNameMap");
		exit(1);
	}
}
/************************************************************/



