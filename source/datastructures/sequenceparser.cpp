/*
 *  sequenceParser.cpp
 *  Mothur
 *
 *  Created by westcott on 9/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "sequenceparser.h"
#include "sharedutilities.h"

/************************************************************/
SequenceParser::SequenceParser(string groupFile, string fastaFile, string nameFile, vector<string> groupsSelected) {
	try {
		
		m = MothurOut::getInstance();
		int error;
		
		//read group file
		groupMap = new GroupMap(groupFile);
		error = groupMap->readMap();
		
		if (error == 1) { m->setControl_pressed(true); }
		
		//initialize maps
        vector<string> namesOfGroups = groupMap->getNamesOfGroups();
        set<string> selectedGroups;
        if (groupsSelected.size() != 0) {
            SharedUtil util;  util.setGroups(groupsSelected, namesOfGroups);
            namesOfGroups = groupsSelected;
        }
        
        for (int i = 0; i < namesOfGroups.size(); i++) {
            vector<Sequence> temp;
            map<string, string> tempMap;
            seqs[namesOfGroups[i]] = temp;
            nameMapPerGroup[namesOfGroups[i]] = tempMap;
            selectedGroups.insert(namesOfGroups[i]);
        }
        
        map<string, string>::iterator it;
        map<string, string> nameMap;
        m->readNames(nameFile, nameMap);
        
		//read fasta file making sure each sequence is in the group file
		ifstream in;
		m->openInputFile(fastaFile, in);
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            Sequence seq(in); m->gobble(in);
            
            if (seq.getName() != "") {
                
                it = nameMap.find(seq.getName());
                
                if (it != nameMap.end()) { //in namefile
                    
                    vector<string> names;
                    string secondCol = it->second;
                    m->splitAtChar(secondCol, names, ',');
                    
                    map<string, string> splitMap; //group -> name1,name2,...
                    map<string, string>::iterator itSplit;
                    for (int i = 0; i < names.size(); i++) {
                        string group = groupMap->getGroup(names[i]);
                        
                        if (selectedGroups.count(group) != 0) { //this is a group we want
                            if (group == "not found") {  error = 1; m->mothurOut("[ERROR]: " + names[i] + " is in your names file and not in your group file, please correct.\n");  }
                            else {
                                allSeqsMap[names[i]] = names[0];
                                
                                itSplit = splitMap.find(group);
                                if (itSplit != splitMap.end()) { //adding seqs to this group
                                    (itSplit->second) += "," + names[i];
                                }else { //first sighting of this group
                                    splitMap[group] = names[i];
                                }
                            }
                        }
                    }
                    
                    //fill nameMapPerGroup - holds all lines in namefile separated by group
                    for (itSplit = splitMap.begin(); itSplit != splitMap.end(); itSplit++) {
                        //grab first name
                        string firstName = "";
                        for(int i = 0; i < (itSplit->second).length(); i++) {
                            if (((itSplit->second)[i]) != ',') {
                                firstName += ((itSplit->second)[i]);
                            }else { break; }
                        }
                        
                        //group1 -> seq1 -> seq1,seq2,seq3
                        nameMapPerGroup[itSplit->first][firstName] = itSplit->second;
                        seqs[itSplit->first].push_back(Sequence(firstName, seq.getAligned()));
                    }

                }else { error = 1; m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your name file, please correct.\n"); }
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
SequenceParser::SequenceParser(string groupFile, string fastaFile, vector<string> groupsSelected) {
	try {
		
		m = MothurOut::getInstance();
		int error;
		
		//read group file
		groupMap = new GroupMap(groupFile);
		error = groupMap->readMap();
		
		if (error == 1) { m->setControl_pressed(true); }
		
		//initialize maps
        vector<string> namesOfGroups = groupMap->getNamesOfGroups();
        set<string> selectedGroups;
        if (groupsSelected.size() != 0) {
            SharedUtil util;  util.setGroups(groupsSelected, namesOfGroups);
            namesOfGroups = groupsSelected;
        }
        
		for (int i = 0; i < namesOfGroups.size(); i++) {
			vector<Sequence> temp;
			seqs[namesOfGroups[i]] = temp;
            selectedGroups.insert(namesOfGroups[i]);
		}
		
		//read fasta file making sure each sequence is in the group file
		ifstream in;
		m->openInputFile(fastaFile, in);
		
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { break; }
			
			Sequence seq(in); m->gobble(in);
			
			if (seq.getName() != "") {
                
				string group = groupMap->getGroup(seq.getName());
                if (selectedGroups.count(group) != 0) { //this is a group we want
                    if (group == "not found") {  error = 1; m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your groupfile, please correct."); m->mothurOutEndLine();  }
                    else {
                        seqs[group].push_back(seq);
                    }
                }
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
SequenceParser::~SequenceParser(){ delete groupMap; }
/************************************************************/
int SequenceParser::getNumGroups(){ return groupMap->getNumGroups(); }
/************************************************************/
vector<string> SequenceParser::getNamesOfGroups(){ return groupMap->getNamesOfGroups(); }
/************************************************************/
bool SequenceParser::isValidGroup(string g){ return groupMap->isValidGroup(g); }
/************************************************************/
int SequenceParser::getNumSeqs(string g){ 
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
		m->errorOut(e, "SequenceParser", "getNumSeqs");
		exit(1);
	}
}
/************************************************************/
vector<Sequence> SequenceParser::getSeqs(string g){ 
	try {
		map<string, vector<Sequence> >::iterator it;
		vector<Sequence> seqForThisGroup;
		
		it = seqs.find(g);
		if(it == seqs.end()) {
			m->mothurOut("[ERROR]: No sequences available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			seqForThisGroup = it->second;
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: group " + g + " fasta file has " + toString(seqForThisGroup.size()) + " sequences.");  }
		}
		
		return seqForThisGroup; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getSeqs");
		exit(1);
	}
}
/************************************************************/
int SequenceParser::getSeqs(string g, string filename, string tag, string tag2, long long& numSeqs, bool uchimeFormat=false){
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
						int num = m->getNumNames(itNameMap->second);
						
						seqPriorityNode temp(num, seqForThisGroup[i].getUnaligned(), seqForThisGroup[i].getName());
						nameVector.push_back(temp);
					}
				}
				
				if (error == 1) { out.close(); m->mothurRemove(filename); return 1; }
				
				//sort by num represented
				sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);

				//print new file in order of
				for (int i = 0; i < nameVector.size(); i++) {
					
					if(m->getControl_pressed()) { out.close(); m->mothurRemove(filename); return 1; }
					
					out << ">" <<  nameVector[i].name << tag << nameVector[i].numIdentical << tag2 << endl << nameVector[i].seq << endl; //
				}
				
			}else { 
                //m->mothurOut("Group " + g +  " contains " + toString(seqForThisGroup.size()) + " unique seqs.\n");
				for (int i = 0; i < seqForThisGroup.size(); i++) {
					
					if(m->getControl_pressed()) { out.close(); m->mothurRemove(filename); return 1; }
					
					seqForThisGroup[i].printSequence(out);	
				}
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
map<string, string> SequenceParser::getNameMap(string g){ 
	try {
		map<string, map<string, string> >::iterator it;
		map<string, string> nameMapForThisGroup;
		
		it = nameMapPerGroup.find(g);
		if(it == nameMapPerGroup.end()) {
			m->mothurOut("[ERROR]: No nameMap available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			nameMapForThisGroup = it->second;
            if (m->getDebug()) {  m->mothurOut("[DEBUG]: group " + g + " name file has " + toString(nameMapForThisGroup.size()) + " unique sequences.");  }
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
		map<string, map<string, string> >::iterator it;
		map<string, string> nameMapForThisGroup;
		
		it = nameMapPerGroup.find(g);
		if(it == nameMapPerGroup.end()) {
			m->mothurOut("[ERROR]: No nameMap available for group " + g + ", please correct."); m->mothurOutEndLine();
		}else {
			nameMapForThisGroup = it->second;
			
			ofstream out;
			m->openOutputFile(filename, out);
			
			for (map<string, string>::iterator itFile = nameMapForThisGroup.begin(); itFile != nameMapForThisGroup.end(); itFile++) {
				
				if(m->getControl_pressed()) { out.close(); m->mothurRemove(filename); return 1; }
				
				out << itFile->first << '\t' << itFile->second << endl;
			}
			
			out.close();
		}
		
		return 0; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getNameMap");
		exit(1);
	}
}
/************************************************************/



