/*
 *  sequenceParser.cpp
 *  Mothur
 *
 *  Created by westcott on 9/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "sequenceParser.h"


/************************************************************/
SequenceParser::SequenceParser(string groupFile, string fastaFile, string nameFile) {
	try {
		
		m = MothurOut::getInstance();
		int error;
		
		//read group file
		groupMap = new GroupMap(groupFile);
		error = groupMap->readMap();
		
		if (error == 1) { m->control_pressed = true; }
		
		//initialize maps
		vector<string> namesOfGroups = groupMap->getNamesOfGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {
			vector<Sequence> temp;
			map<string, string> tempMap;
			seqs[namesOfGroups[i]] = temp;
			nameMapPerGroup[namesOfGroups[i]] = tempMap;
		}
		
		//read fasta file making sure each sequence is in the group file
		ifstream in;
		m->openInputFile(fastaFile, in);
		
		map<string, string> seqName; //stores name -> sequence string so we can make new "unique" sequences when we parse the name file
		while (!in.eof()) {
			
			if (m->control_pressed) { break; }
			
			Sequence seq(in); m->gobble(in);
			
			if (seq.getName() != "") {
				
				 string group = groupMap->getGroup(seq.getName());
				 if (group == "not found") {  error = 1; m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your groupfile, please correct."); m->mothurOutEndLine();  }
				 else {	
					 seqs[group].push_back(seq);	
					 seqName[seq.getName()] = seq.getAligned();
				 }
			}
		}
		in.close();
				 
		if (error == 1) { m->control_pressed = true; }
				 
		//read name file
		ifstream inName;
		m->openInputFile(nameFile, inName);
		
		string first, second;
		int countName = 0;
		while(!inName.eof()) {
			
			if (m->control_pressed) { break; }
			
			inName >> first; m->gobble(inName);
			inName >> second; m->gobble(inName);
			
			vector<string> names;
			m->splitAtChar(second, names, ',');
			
			//get aligned string for these seqs from the fasta file
			string alignedString = "";
			map<string, string>::iterator itAligned = seqName.find(names[0]);
			if (itAligned == seqName.end()) {
				error = 1; m->mothurOut("[ERROR]: " + names[0] + " is in your name file and not in your fasta file, please correct."); m->mothurOutEndLine();
			}else {
				alignedString = itAligned->second;
			}
			
			//separate by group - parse one line in name file
			map<string, string> splitMap; //group -> name1,name2,...
			map<string, string>::iterator it;
			for (int i = 0; i < names.size(); i++) {
				
				string group = groupMap->getGroup(names[i]);
				if (group == "not found") {  error = 1; m->mothurOut("[ERROR]: " + names[i] + " is in your name file and not in your groupfile, please correct."); m->mothurOutEndLine();  }
				else {	
					
					it = splitMap.find(group);
					if (it != splitMap.end()) { //adding seqs to this group
						(it->second) += "," + names[i];
						countName++;
					}else { //first sighting of this group
						splitMap[group] = names[i];
						countName++;
						
						//is this seq in the fasta file?
						if (i != 0) { //if not then we need to add a duplicate sequence to the seqs for this group so the new "fasta" and "name" files will match
							Sequence tempSeq(names[i], alignedString); //get the first guys sequence string since he's in the fasta file.
							seqs[group].push_back(tempSeq);
						}
					}
				}
				
				allSeqsMap[names[i]] = names[0];
			}
			
			
			//fill nameMapPerGroup - holds all lines in namefile separated by group
			for (it = splitMap.begin(); it != splitMap.end(); it++) {
				//grab first name
				string firstName = "";
				for(int i = 0; i < (it->second).length(); i++) {
					if (((it->second)[i]) != ',') {
						firstName += ((it->second)[i]);
					}else { break; }
				}
				
				//group1 -> seq1 -> seq1,seq2,seq3
				nameMapPerGroup[it->first][firstName] = it->second;
			}
		}
		
		inName.close();
		
		if (error == 1) { m->control_pressed = true; }
		
		if (countName != (groupMap->getNumSeqs())) {
			m->mothurOutEndLine();
			m->mothurOut("[ERROR]: Your name file contains " + toString(countName) + " valid sequences, and your groupfile contains " + toString(groupMap->getNumSeqs()) + ", please correct.");
			m->mothurOutEndLine();
			m->control_pressed = true;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "SequenceParser");
		exit(1);
	}
}
/************************************************************/
SequenceParser::SequenceParser(string groupFile, string fastaFile) {
	try {
		
		m = MothurOut::getInstance();
		int error;
		
		//read group file
		groupMap = new GroupMap(groupFile);
		error = groupMap->readMap();
		
		if (error == 1) { m->control_pressed = true; }
		
		//initialize maps
		vector<string> namesOfGroups = groupMap->getNamesOfGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {
			vector<Sequence> temp;
			seqs[namesOfGroups[i]] = temp;
		}
		
		//read fasta file making sure each sequence is in the group file
		ifstream in;
		m->openInputFile(fastaFile, in);
		int count = 0;
		
		while (!in.eof()) {
			
			if (m->control_pressed) { break; }
			
			Sequence seq(in); m->gobble(in);
			
			if (seq.getName() != "") {
				
				string group = groupMap->getGroup(seq.getName());
				if (group == "not found") {  error = 1; m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file and not in your groupfile, please correct."); m->mothurOutEndLine();  }
				else {	seqs[group].push_back(seq);	count++; }
			}
		}
		in.close();
		
		if (error == 1) { m->control_pressed = true; }
		
		if (count != (groupMap->getNumSeqs())) {
			m->mothurOutEndLine();
			m->mothurOut("[ERROR]: Your fasta file contains " + toString(count) + " valid sequences, and your groupfile contains " + toString(groupMap->getNumSeqs()) + ", please correct.");
			if (count < (groupMap->getNumSeqs())) { m->mothurOut(" Did you forget to include the name file?"); }
			m->mothurOutEndLine();
			m->control_pressed = true;
		}
		
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
string SequenceParser::getGroup(string g){ return groupMap->getGroup(g); }
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
		}
		
		return seqForThisGroup; 
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceParser", "getSeqs");
		exit(1);
	}
}
/************************************************************/
int SequenceParser::getSeqs(string g, string filename, bool uchimeFormat=false){ 
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
						
						seqPriorityNode temp(num, seqForThisGroup[i].getAligned(), seqForThisGroup[i].getName());
						nameVector.push_back(temp);
					}
				}
				
				if (error == 1) { out.close(); m->mothurRemove(filename); return 1; }
				
				//sort by num represented
				sort(nameVector.begin(), nameVector.end(), compareSeqPriorityNodes);

				//print new file in order of
				for (int i = 0; i < nameVector.size(); i++) {
					
					if(m->control_pressed) { out.close(); m->mothurRemove(filename); return 1; }
					
					out << ">" << nameVector[i].name  << "/ab=" << nameVector[i].numIdentical << "/" << endl << nameVector[i].seq << endl;
				}
				
			}else { 
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
				
				if(m->control_pressed) { out.close(); m->mothurRemove(filename); return 1; }
				
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



