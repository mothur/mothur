/*
 *  splitmatrix.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitmatrix.h"
#include "phylotree.h"

/***********************************************************************/

SplitMatrix::SplitMatrix(string distfile, string name, string tax, float c, string t, bool l){
	m = MothurOut::getInstance();
	distFile = distfile;
	cutoff = c;
	namefile = name;
	method = t;
	taxFile = tax;
	large = l;
}

/***********************************************************************/

int SplitMatrix::split(){
	try {
        
		if (method == "distance") {  
			splitDistance();
		}else if (method == "classify") {
			splitClassify();
		}else {
			m->mothurOut("Unknown splitting method, aborting split."); m->mothurOutEndLine();
			map<string, string> temp;
			temp[distFile] = namefile;
			dists.push_back(temp);
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "split");
		exit(1);
	}
}
/***********************************************************************/
int SplitMatrix::splitDistance(){
	try {
        
		if (large)	{ splitDistanceLarge(); }
		else		{ splitDistanceRAM();	}
			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitDistance");
		exit(1);
	}
}

/***********************************************************************/
int SplitMatrix::splitClassify(){
	try {
		cutoff = int(cutoff);
		
		map<string, int> seqGroup;
		map<string, int>::iterator it;
		map<string, int>::iterator it2;
		
		int numGroups = 0;
		
		//build tree from users taxonomy file
		PhyloTree* phylo = new PhyloTree();
		
		ifstream in;
		openInputFile(taxFile, in);
			
		//read in users taxonomy file and add sequences to tree
		string seqname, tax;
		while(!in.eof()){
			in >> seqname >> tax; gobble(in);
				
			phylo->addSeqToTree(seqname, tax);
		}
		in.close();
		
		phylo->assignHeirarchyIDs(0);

		//make sure the cutoff is not greater than maxlevel
		if (cutoff > phylo->getMaxLevel()) { m->mothurOut("splitcutoff is greater than the longest taxonomy, using " + toString(phylo->getMaxLevel())); m->mothurOutEndLine(); cutoff = phylo->getMaxLevel(); }
		
		//for each node in tree
		for (int i = 0; i < phylo->getNumNodes(); i++) {
		
			//is this node within the cutoff
			TaxNode taxon = phylo->get(i);
		
			if (taxon.level == cutoff) {//if yes, then create group containing this nodes sequences
				if (taxon.accessions.size() > 1) { //if this taxon just has one seq its a singleton
					for (int j = 0; j < taxon.accessions.size(); j++) {
						seqGroup[taxon.accessions[j]] = numGroups;
					}
					numGroups++;
				}
			}
		}

		ifstream dFile;
		openInputFile(distFile, dFile);
		ofstream outFile;
		
		for (int i = 0; i < numGroups; i++) { //remove old temp files, just in case
			remove((distFile + "." + toString(i) + ".temp").c_str());
		}
		
		
		//for buffering the io to improve speed
		 //allow for 10 dists to be stored, then output.
		vector<string> outputs;  outputs.resize(numGroups, "");
		vector<int> numOutputs;	 numOutputs.resize(numGroups, 0);	
		
		//for each distance
		while(dFile){
			string seqA, seqB;
			float dist;
			
			if (m->control_pressed) { dFile.close(); for (int i = 0; i < numGroups; i++) { remove((distFile + "." + toString(i) + ".temp").c_str());	} }
			
			dFile >> seqA >> seqB >> dist;  gobble(dFile);
			
			//if both sequences are in the same group then they are within the cutoff
			it = seqGroup.find(seqA);
			it2 = seqGroup.find(seqB);
			
			if ((it != seqGroup.end()) && (it2 != seqGroup.end())) { //they are both not singletons 
				if (it->second == it2->second) { //they are from the same group so add the distance
					if (numOutputs[it->second] > 10) {
						openOutputFileAppend((distFile + "." + toString(it->second) + ".temp"), outFile);
						outFile << outputs[it->second] << seqA << '\t' << seqB << '\t' << dist << endl;
						outFile.close();
						outputs[it->second] = "";
						numOutputs[it->second] = 0;
					}else{
						outputs[it->second] += seqA + '\t' + seqB + '\t' + toString(dist)  + '\n';
						numOutputs[it->second]++;
					}
				}
			}
		}
		dFile.close();
	
		for (int i = 0; i < numGroups; i++) { //remove old temp files, just in case
			remove((namefile + "." + toString(i) + ".temp").c_str());
			
			//write out any remaining buffers
			if (numOutputs[it->second] > 0) {
				openOutputFileAppend((distFile + "." + toString(i) + ".temp"), outFile);
				outFile << outputs[i];
				outFile.close();
				outputs[i] = "";
				numOutputs[i] = 0;
			}
		}
		
		ifstream bigNameFile;
		openInputFile(namefile, bigNameFile);
		
		singleton = namefile + ".extra.temp";
		ofstream remainingNames;
		openOutputFile(singleton, remainingNames);
		
		bool wroteExtra = false;
						
		string name, nameList;
		while(!bigNameFile.eof()){
			bigNameFile >> name >> nameList;  gobble(bigNameFile);
			
			//did this sequence get assigned a group
			it = seqGroup.find(name);
			
			if (it != seqGroup.end()) {  
				openOutputFileAppend((namefile + "." + toString(it->second) + ".temp"), outFile);
				outFile << name << '\t' << nameList << endl;
				outFile.close();
			}else{
				wroteExtra = true;
				remainingNames << name << '\t' << nameList << endl;
			}
		}
		bigNameFile.close();
		remainingNames.close();
		
		if (!wroteExtra) { 
			remove(singleton.c_str());
			singleton = "none";
		}
			
		for(int i=0;i<numGroups;i++){
			string tempNameFile = namefile + "." + toString(i) + ".temp";
			string tempDistFile = distFile + "." + toString(i) + ".temp";
				
			map<string, string> temp;
			temp[tempDistFile] = tempNameFile;
			dists.push_back(temp);
		}
		
		if (m->control_pressed)	 {  
			for (int i = 0; i < dists.size(); i++) { 
				remove((dists[i].begin()->first).c_str());
				remove((dists[i].begin()->second).c_str());
			}
			dists.clear();
		}
		
		return 0;
			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitClassify");
		exit(1);
	}
}
/***********************************************************************/
int SplitMatrix::splitDistanceLarge(){
	try {
		vector<set<string> > groups;
		
		//for buffering the io to improve speed
		 //allow for 30 dists to be stored, then output.
		vector<string> outputs;
		vector<int> numOutputs;
		vector<bool> wroteOutPut;
		
		int numGroups = 0;

		ofstream outFile;
		ifstream dFile;
		openInputFile(distFile, dFile);
	
		while(dFile){
			string seqA, seqB;
			float dist;

			dFile >> seqA >> seqB >> dist;
			
			if (m->control_pressed) {   dFile.close();  for(int i=0;i<numGroups;i++){	if(groups[i].size() > 0){  remove((distFile + "." + toString(i) + ".temp").c_str()); }  } return 0; }
					
			if(dist < cutoff){
				//cout << "in cutoff: " << dist << endl;
				int groupIDA = -1;
				int groupIDB = -1;
				int groupID = -1;
				
				for(int i=0;i<numGroups;i++){
					set<string>::iterator aIt = groups[i].find(seqA);
					set<string>::iterator bIt = groups[i].find(seqB);
					
					if(groupIDA == -1 && aIt != groups[i].end()){//seqA is not already assigned to a group and is in group[i], so assign seqB to group[i]
						groups[i].insert(seqB);
						groupIDA = i;
						groupID = groupIDA;

						//cout << "in aIt: " << groupID << endl;
	//					break;
					}
					else if(groupIDB == -1 && bIt != groups[i].end()){//seqB is not already assigned to a group and is in group[i], so assign seqA to group[i]
						groups[i].insert(seqA);
						groupIDB = i;
						groupID = groupIDB;

					//	cout << "in bIt: " << groupID << endl;
	//					break;
					}
				
					if(groupIDA != -1 && groupIDB != -1){//both ifs above have been executed, so we need to decide who to assign them to
						if(groupIDA < groupIDB){
						//	cout << "A: " << groupIDA << "\t" << groupIDB << endl;
							groups[groupIDA].insert(groups[groupIDB].begin(), groups[groupIDB].end()); //merge two groups into groupIDA
							groups[groupIDB].clear(); 
							groupID = groupIDA;
						}
						else{
						//	cout << "B: " << groupIDA << "\t" << groupIDB << endl;
							groups[groupIDB].insert(groups[groupIDA].begin(), groups[groupIDA].end()); //merge two groups into groupIDB
							groups[groupIDA].clear();  
							groupID = groupIDB;
						}
						break;
					}
				}
				
	//windows is gonna gag on the reuse of outFile, will need to make it local...
				
				if(groupIDA == -1 && groupIDB == -1){ //we need a new group
					set<string> newGroup;
					newGroup.insert(seqA);
					newGroup.insert(seqB);
					groups.push_back(newGroup);
									
					string tempOut = seqA + '\t' + seqB + '\t' + toString(dist) + '\n';
					outputs.push_back(tempOut);
					numOutputs.push_back(1);
					wroteOutPut.push_back(false);
					
					numGroups++;
				}
				else{
					string fileName = distFile + "." + toString(groupID) + ".temp";
											
					//have we reached the max buffer size
					if (numOutputs[groupID] > 60) { //write out sequence
						outFile.open(fileName.c_str(), ios::app);
						outFile << outputs[groupID] << seqA << '\t' << seqB << '\t' << dist << endl;
						outFile.close();
						
						outputs[groupID] = "";
						numOutputs[groupID] = 0;
						wroteOutPut[groupID] = true;
					}else {
						outputs[groupID] +=  seqA + '\t' + seqB + '\t' + toString(dist)  + '\n';
						numOutputs[groupID]++;
					}
					
					if(groupIDA != -1 && groupIDB != -1){ //merge distance files of two groups you merged above
						string row, column, distance;
						if(groupIDA<groupIDB){
							
							//merge memory
							numOutputs[groupID] += numOutputs[groupIDB];
							outputs[groupID] += outputs[groupIDB];
							
							outputs[groupIDB] = "";
							numOutputs[groupIDB] = 0;
							
							//if groupB is written to file it is above buffer size so read and write to new merged file
							if (wroteOutPut[groupIDB]) {
								string fileName2 = distFile + "." + toString(groupIDB) + ".temp";
								ifstream fileB(fileName2.c_str(), ios::ate);
								
								outFile.open(fileName.c_str(), ios::app);
								
								long size;
								char* memblock;

								size = fileB.tellg();
				
								fileB.seekg (0, ios::beg);
								
								int numRead = size / 1024;
								int lastRead = size % 1024;

								for (int i = 0; i < numRead; i++) {
				
									memblock = new char [1024];
								
									fileB.read (memblock, 1024);
									
									string temp = memblock;
									outFile << temp.substr(0, 1024);
									
									delete memblock;
								}
								
								memblock = new char [lastRead];
								
								fileB.read (memblock, lastRead);
								
								//not sure why but it will read more than lastRead char...??
								string temp = memblock;
								outFile << temp.substr(0, lastRead);
								delete memblock;
								
								fileB.close();
								remove(fileName2.c_str());
								
								//write out the merged memory
								if (numOutputs[groupID] > 60) {
									outFile << outputs[groupID];
									outputs[groupID] = "";
									numOutputs[groupID] = 0;
								}
								
								outFile.close();
								
								wroteOutPut[groupID] = true;
								wroteOutPut[groupIDB] = false;
							}else{ } //just merge b's memory with a's memory 
						}
						else{
							numOutputs[groupID] += numOutputs[groupIDA];
							outputs[groupID] += outputs[groupIDA];
							
							outputs[groupIDA] = "";
							numOutputs[groupIDA] = 0;
							
							if (wroteOutPut[groupIDA]) {
								string fileName2 = distFile + "." + toString(groupIDA) + ".temp";
								ifstream fileB(fileName2.c_str(), ios::ate);
								
								outFile.open(fileName.c_str(), ios::app);
								
								long size;
								char* memblock;

								size = fileB.tellg();
															
								fileB.seekg (0, ios::beg);
								
								int numRead = size / 1024;
								int lastRead = size % 1024;

								for (int i = 0; i < numRead; i++) {
				
									memblock = new char [1024];
								
									fileB.read (memblock, 1024);
									string temp = memblock;
									outFile << temp.substr(0, 1024);
									
									delete memblock;
								}
								
								memblock = new char [lastRead];
								
								fileB.read (memblock, lastRead);
								
								//not sure why but it will read more than lastRead char...??
								string temp = memblock;
								outFile << temp.substr(0, lastRead);
									
								delete memblock;
								
								fileB.close();
								remove(fileName2.c_str());
								
								//write out the merged memory
								if (numOutputs[groupID] > 60) {
									outFile << outputs[groupID];
									outputs[groupID] = "";
									numOutputs[groupID] = 0;
								}
								
								outFile.close();
								
								wroteOutPut[groupID] = true;
								wroteOutPut[groupIDA] = false;
							}else { } //just merge memory
						}					
					}
				}
			}
			gobble(dFile);
		}
		dFile.close();
		
		for (int i = 0; i < numGroups; i++) {
			if (numOutputs[i] > 0) {
				string fileName = distFile + "." + toString(i) + ".temp";
				outFile.open(fileName.c_str(), ios::app);
				outFile << outputs[i];
				outFile.close();
			}
		}

		splitNames(groups);
				
		return 0;			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitDistanceLarge");
		exit(1);
	}
}
//********************************************************************************************************************
int SplitMatrix::splitNames(vector<set<string> >& groups){
	try {
		int numGroups = groups.size();
	
		ifstream bigNameFile(namefile.c_str());
		if(!bigNameFile){
			cerr << "Error: We can't open the name file\n";
			exit(1);
		}
		
		map<string, string> nameMap;
		string name, nameList;
		while(bigNameFile){
			bigNameFile >> name >> nameList;
			nameMap[name] = nameList;
			gobble(bigNameFile);
		}
		bigNameFile.close();
			
		for(int i=0;i<numGroups;i++){  //parse names file to match distance files
			int numSeqsInGroup = groups[i].size();
			
			if(numSeqsInGroup > 0){
				string fileName = namefile + "." + toString(i) + ".temp";
				ofstream smallNameFile(fileName.c_str(), ios::ate);
				
				for(set<string>::iterator gIt=groups[i].begin();gIt!=groups[i].end();gIt++){
					map<string,string>::iterator nIt = nameMap.find(*gIt);
					if (nIt != nameMap.end()) {
						smallNameFile << nIt->first << '\t' << nIt->second << endl;
						nameMap.erase(nIt);
					}else{
						m->mothurOut((*gIt) + " is in your distance file and not in your namefile.  Please correct."); m->mothurOutEndLine(); exit(1);
					}
				}
				smallNameFile.close();
			}
		}
		
		//names of singletons
		if (nameMap.size() != 0) {
			singleton = namefile + ".extra.temp";
			ofstream remainingNames(singleton.c_str(), ios::ate);
			for(map<string,string>::iterator nIt=nameMap.begin();nIt!=nameMap.end();nIt++){
				remainingNames << nIt->first << '\t' << nIt->second << endl;
			}
			remainingNames.close();
		}else { singleton = "none"; }
			
		for(int i=0;i<numGroups;i++){
			if(groups[i].size() > 0){
				string tempNameFile = namefile + "." + toString(i) + ".temp";
				string tempDistFile = distFile + "." + toString(i) + ".temp";
				
				map<string, string> temp;
				temp[tempDistFile] = tempNameFile;
				dists.push_back(temp);
			}
		}
		
		if (m->control_pressed)	 {  
			for (int i = 0; i < dists.size(); i++) { 
				remove((dists[i].begin()->first).c_str());
				remove((dists[i].begin()->second).c_str());
			}
			dists.clear();
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitNames");
		exit(1);
	}
}
//********************************************************************************************************************
int SplitMatrix::splitDistanceRAM(){
	try {
		vector<set<string> > groups;
		vector<string> outputs;
		
		int numGroups = 0;

		ifstream dFile;
		openInputFile(distFile, dFile);

		while(dFile){
			string seqA, seqB;
			float dist;

			dFile >> seqA >> seqB >> dist;
			
			if (m->control_pressed) {   dFile.close();  for(int i=0;i<numGroups;i++){	if(groups[i].size() > 0){  remove((distFile + "." + toString(i) + ".temp").c_str()); }  } return 0; }
					
			if(dist < cutoff){
				//cout << "in cutoff: " << dist << endl;
				int groupIDA = -1;
				int groupIDB = -1;
				int groupID = -1;
				
				for(int i=0;i<numGroups;i++){
					set<string>::iterator aIt = groups[i].find(seqA);
					set<string>::iterator bIt = groups[i].find(seqB);
					
					if(groupIDA == -1 && aIt != groups[i].end()){//seqA is not already assigned to a group and is in group[i], so assign seqB to group[i]
						groups[i].insert(seqB);
						groupIDA = i;
						groupID = groupIDA;

						//cout << "in aIt: " << groupID << endl;
	//					break;
					}
					else if(groupIDB == -1 && bIt != groups[i].end()){//seqB is not already assigned to a group and is in group[i], so assign seqA to group[i]
						groups[i].insert(seqA);
						groupIDB = i;
						groupID = groupIDB;

					//	cout << "in bIt: " << groupID << endl;
	//					break;
					}
				
					if(groupIDA != -1 && groupIDB != -1){//both ifs above have been executed, so we need to decide who to assign them to
						if(groupIDA < groupIDB){
						//	cout << "A: " << groupIDA << "\t" << groupIDB << endl;
							groups[groupIDA].insert(groups[groupIDB].begin(), groups[groupIDB].end()); //merge two groups into groupIDA
							groups[groupIDB].clear(); 
							groupID = groupIDA;
						}
						else{
						//	cout << "B: " << groupIDA << "\t" << groupIDB << endl;
							groups[groupIDB].insert(groups[groupIDA].begin(), groups[groupIDA].end()); //merge two groups into groupIDB
							groups[groupIDA].clear();  
							groupID = groupIDB;
						}
						break;
					}
				}
				
	//windows is gonna gag on the reuse of outFile, will need to make it local...
				
				if(groupIDA == -1 && groupIDB == -1){ //we need a new group
					set<string> newGroup;
					newGroup.insert(seqA);
					newGroup.insert(seqB);
					groups.push_back(newGroup);
									
					string tempOut = seqA + '\t' + seqB + '\t' + toString(dist) + '\n';
					outputs.push_back(tempOut);
					numGroups++;
				}
				else{
											
					outputs[groupID] +=  seqA + '\t' + seqB + '\t' + toString(dist)  + '\n';
					
					if(groupIDA != -1 && groupIDB != -1){ //merge distance files of two groups you merged above
						string row, column, distance;
						if(groupIDA<groupIDB){
							//merge memory
							outputs[groupID] += outputs[groupIDB];
							outputs[groupIDB] = "";
						}else{
							outputs[groupID] += outputs[groupIDA];
							outputs[groupIDA] = "";
						}					
					}
				}
			}
			gobble(dFile);
		}
		dFile.close();
		
		for (int i = 0; i < numGroups; i++) {
			if (outputs[i] != "") {
				ofstream outFile;
				string fileName = distFile + "." + toString(i) + ".temp";
				outFile.open(fileName.c_str(), ios::ate);
				outFile << outputs[i];
				outFile.close();
			}
		}

		splitNames(groups);
				
		return 0;			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitDistanceRAM");
		exit(1);
	}
}
//********************************************************************************************************************
//sorts biggest to smallest
inline bool compareFileSizes(map<string, string> left, map<string, string> right){
	
	FILE * pFile;
	long leftsize = 0;
		
	//get num bytes in file
	string filename = left.begin()->first;
	pFile = fopen (filename.c_str(),"rb");
	string error = "Error opening " + filename;
	if (pFile==NULL) perror (error.c_str());
	else{
		fseek (pFile, 0, SEEK_END);
		leftsize=ftell (pFile);
		fclose (pFile);
	}

	FILE * pFile2;
	long rightsize = 0;
		
	//get num bytes in file
	filename = right.begin()->first;
	pFile2 = fopen (filename.c_str(),"rb");
	error = "Error opening " + filename;
	if (pFile2==NULL) perror (error.c_str());
	else{
		fseek (pFile2, 0, SEEK_END);
		rightsize=ftell (pFile2);
		fclose (pFile2);
	}

	return (leftsize > rightsize);	
} 
/***********************************************************************/
//returns map of distance files -> namefile sorted by distance file size
vector< map< string, string> > SplitMatrix::getDistanceFiles(){
	try {	
		
		sort(dists.begin(), dists.end(), compareFileSizes);
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "getDistanceFiles");
		exit(1);
	}
}
/***********************************************************************/
SplitMatrix::~SplitMatrix(){}
/***********************************************************************/

