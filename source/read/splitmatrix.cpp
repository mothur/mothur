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
#include "distancecommand.h"
#include "seqsummarycommand.h"

/***********************************************************************/

SplitMatrix::SplitMatrix(string distfile, string name, string count, string tax, float c, string t, bool l){
	m = MothurOut::getInstance();
	distFile = distfile;
	cutoff = c;
	namefile = name;
	method = t;
	taxFile = tax;
    countfile = count;
	large = l;
}
/***********************************************************************/

SplitMatrix::SplitMatrix(string ffile, string name, string count, string tax, float c, float cu, string t, int p, bool cl, string output){
	m = MothurOut::getInstance();
	fastafile = ffile;
	namefile = name;
    countfile = count;
	taxFile = tax;
	cutoff = c;  //tax level cutoff
	distCutoff = cu; //for fasta method if you are creating distance matrix you need a cutoff for that
	method = t;
	processors = p;
    classic = cl;
	outputDir = output;
}

/***********************************************************************/

int SplitMatrix::split(){
	try {
        
		if (method == "distance") {  
			splitDistance();
		}else if ((method == "classify") || (method == "fasta")) {
			splitClassify();
		}else {
			m->mothurOut("Unknown splitting method, aborting split."); m->mothurOutEndLine();
			map<string, string> temp;
			if (namefile != "") {  temp[distFile] = namefile; }
            else { temp[distFile] = countfile; }
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
		
		return 0;
			
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
		
        map<string, string> temp;
        m->readTax(taxFile, temp, true);
        
        for (map<string, string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {
            phylo->addSeqToTree(itTemp->first, itTemp->second);
            temp.erase(itTemp++);
        }
		
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
	
		delete phylo;
		
		if (method == "classify") {
			splitDistanceFileByTax(seqGroup, numGroups);
		}else {
			createDistanceFilesFromTax(seqGroup, numGroups);
		}
		
		return 0;
			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitClassify");
		exit(1);
	}
}
/***********************************************************************/
int SplitMatrix::createDistanceFilesFromTax(map<string, int>& seqGroup, int numGroups){
	try {
        
		map<string, int> copyGroups = seqGroup;
		map<string, int>::iterator it;
		set<string> names;
				
		ifstream in;
		m->openInputFile(fastafile, in);

		for (int i = 0; i < numGroups; i++) {  m->mothurRemove((fastafile + "." + toString(i) + ".temp")); }
	
		//parse fastafile
		while (!in.eof()) {
			Sequence query(in); m->gobble(in);
			if (query.getName() != "") {
		
				it = seqGroup.find(query.getName());
				
				//save names in case no namefile is given
				if ((namefile == "") && (countfile == "")) {  names.insert(query.getName()); }
			
				if (it != seqGroup.end()) { //not singleton
                    ofstream outFile;
                    m->openOutputFileAppend((fastafile + "." + toString(it->second) + ".temp"), outFile);
                    query.printSequence(outFile);
                    outFile.close();
					copyGroups.erase(query.getName());
				}
			}
		}
		in.close();
        
        bool error = false;
		//warn about sequence in groups that are not in fasta file
		for(it = copyGroups.begin(); it != copyGroups.end(); it++) {
			m->mothurOut("ERROR: " + it->first + " is missing from your fastafile. This could happen if your taxonomy file is not unique and your fastafile is, or it could indicate and error."); m->mothurOutEndLine();
            error = true;
		}
		copyGroups.clear();
        
        if (error) { exit(1); }
        
		//process each distance file
		for (int i = 0; i < numGroups; i++) { 
			
			string options = "";
            if (classic) { options = "fasta=" + (fastafile + "." + toString(i) + ".temp") + ", processors=" + toString(processors) + ", output=lt"; }
            else { options = "fasta=" + (fastafile + "." + toString(i) + ".temp") + ", processors=" + toString(processors) + ", cutoff=" + toString(distCutoff); }
			if (outputDir != "") { options += ", outputdir=" + outputDir; }
            
            m->mothurCalling = true;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            m->mothurOut("Running command: dist.seqs(" + options + ")"); m->mothurOutEndLine();
            m->mothurCalling = true;

			Command* command = new DistanceCommand(options);
			
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            
			command->execute();
			delete command;
            m->mothurCalling = false;

			m->mothurRemove((fastafile + "." + toString(i) + ".temp"));
			
			//remove old names files just in case
			if (namefile != "") { m->mothurRemove((namefile + "." + toString(i) + ".temp")); }
            else { m->mothurRemove((countfile + "." + toString(i) + ".temp")); }
		}
        
        //restore old fasta file name since dist.seqs overwrites it with the temp files
        m->setFastaFile(fastafile);
        
        vector<string> tempDistFiles;    
        for(int i=0;i<numGroups;i++){
            if (outputDir == "") { outputDir = m->hasPath(fastafile); }
            string tempDistFile = "";
            if (classic) { tempDistFile =  outputDir + m->getRootName(m->getSimpleName((fastafile + "." + toString(i) + ".temp"))) + "phylip.dist";}
            else { tempDistFile = outputDir + m->getRootName(m->getSimpleName((fastafile + "." + toString(i) + ".temp"))) + "dist"; }
            tempDistFiles.push_back(tempDistFile);
        }
        
        splitNames(seqGroup, numGroups, tempDistFiles);
        
		if (m->control_pressed)	 {  for (int i = 0; i < dists.size(); i++) { m->mothurRemove((dists[i].begin()->first)); m->mothurRemove((dists[i].begin()->second)); } dists.clear(); }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
		exit(1);
	}
}
/***********************************************************************/
int SplitMatrix::splitDistanceFileByTax(map<string, int>& seqGroup, int numGroups){
	try {
		map<string, int>::iterator it;
		map<string, int>::iterator it2;
		
        ofstream outFile;
		ifstream dFile;
		m->openInputFile(distFile, dFile);
		
		
		for (int i = 0; i < numGroups; i++) { //remove old temp files, just in case
			m->mothurRemove((distFile + "." + toString(i) + ".temp"));
		}
		
		//for buffering the io to improve speed
		 //allow for 10 dists to be stored, then output.
		vector<string> outputs;  outputs.resize(numGroups, "");
		vector<int> numOutputs;	 numOutputs.resize(numGroups, 0);	
		
		//you can have a group made, but their may be no distances in the file for this group if the taxonomy file and distance file don't match
		//this can occur if we have converted the phylip to column, since we reduce the size at that step by using the cutoff value
		vector<bool> validDistances;   validDistances.resize(numGroups, false); 
		
		//for each distance
		while(dFile){
			string seqA, seqB;
			float dist;
			
			if (m->control_pressed) { dFile.close(); for (int i = 0; i < numGroups; i++) { m->mothurRemove((distFile + "." + toString(i) + ".temp"));	} }
			
			dFile >> seqA >> seqB >> dist;  m->gobble(dFile);
			
			//if both sequences are in the same group then they are within the cutoff
			it = seqGroup.find(seqA);
			it2 = seqGroup.find(seqB);
			
			if ((it != seqGroup.end()) && (it2 != seqGroup.end())) { //they are both not singletons 
				if (it->second == it2->second) { //they are from the same group so add the distance
					if (numOutputs[it->second] > 30) {
						m->openOutputFileAppend((distFile + "." + toString(it->second) + ".temp"), outFile);
						outFile << outputs[it->second] << seqA << '\t' << seqB << '\t' << dist << endl;
						outFile.close();
						outputs[it->second] = "";
						numOutputs[it->second] = 0;
						validDistances[it->second] = true;
					}else{
						outputs[it->second] += seqA + '\t' + seqB + '\t' + toString(dist)  + '\n';
						numOutputs[it->second]++;
					}
				}
			}
		}
		dFile.close();
        
        string inputFile = namefile;
        if (countfile != "") { inputFile = countfile; }
        
        vector<string> tempDistFiles;
		for (int i = 0; i < numGroups; i++) { //remove old temp files, just in case
            string tempDistFile = distFile + "." + toString(i) + ".temp";
            tempDistFiles.push_back(tempDistFile);
			m->mothurRemove((inputFile + "." + toString(i) + ".temp"));
			
			//write out any remaining buffers
			if (numOutputs[i] > 0) {
				m->openOutputFileAppend((distFile + "." + toString(i) + ".temp"), outFile);
				outFile << outputs[i];
				outFile.close();
				outputs[i] = "";
				numOutputs[i] = 0;
				validDistances[i] = true;
			}
		}
		
        splitNames(seqGroup, numGroups, tempDistFiles);
        
		if (m->control_pressed)	 {  
			for (int i = 0; i < dists.size(); i++) { 
				m->mothurRemove((dists[i].begin()->first));
				m->mothurRemove((dists[i].begin()->second));
			}
			dists.clear();
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitDistanceFileByTax");
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

		//ofstream outFile;
		ifstream dFile;
		m->openInputFile(distFile, dFile);
	
		while(dFile){
			string seqA, seqB;
			float dist;

			dFile >> seqA >> seqB >> dist;
			
			if (m->control_pressed) {   dFile.close();  for(int i=0;i<numGroups;i++){	if(groups[i].size() > 0){  m->mothurRemove((distFile + "." + toString(i) + ".temp")); }  } return 0; }
					
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
                        ofstream outFile;
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
                                m->appendFiles(fileName2, fileName);
								m->mothurRemove(fileName2);
                        
								
								//write out the merged memory
								if (numOutputs[groupID] > 60) {
                                    ofstream tempOut;
                                    m->openOutputFile(fileName, tempOut);
									tempOut << outputs[groupID];
									outputs[groupID] = "";
									numOutputs[groupID] = 0;
                                    tempOut.close();
								}
								
								//outFile.close();
								
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
                                m->appendFiles(fileName2, fileName);
								m->mothurRemove(fileName2);
								
								//write out the merged memory
								if (numOutputs[groupID] > 60) {
                                    ofstream tempOut;
                                    m->openOutputFile(fileName, tempOut);
									tempOut << outputs[groupID];
									outputs[groupID] = "";
									numOutputs[groupID] = 0;
                                    tempOut.close();
								}
								
								//outFile.close();
								
								wroteOutPut[groupID] = true;
								wroteOutPut[groupIDA] = false;
							}else { } //just merge memory
						}					
					}
				}
			}
			m->gobble(dFile);
		}
		dFile.close();
        
		vector<string> tempDistFiles;
		for (int i = 0; i < numGroups; i++) {
            string fileName = distFile + "." + toString(i) + ".temp";
            tempDistFiles.push_back(fileName);
            //remove old names files just in case
			
			if (numOutputs[i] > 0) {
                ofstream outFile;
				outFile.open(fileName.c_str(), ios::app);
				outFile << outputs[i];
				outFile.close();
			}
		}
        
        map<string, int> seqGroup;
        for (int i = 0; i < groups.size(); i++) {
            for (set<string>::iterator itNames = groups[i].begin(); itNames != groups[i].end();) {
                seqGroup[*itNames] = i;
                groups[i].erase(itNames++);
            }
        }
        
		splitNames(seqGroup, numGroups, tempDistFiles);
				
		return 0;			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitDistanceLarge");
		exit(1);
	}
}
//********************************************************************************************************************
int SplitMatrix::splitNames(map<string, int>& seqGroup, int numGroups, vector<string>& tempDistFiles){
	try {
        ofstream outFile;
        map<string, int>::iterator it;
        
        string inputFile = namefile;
        if (countfile != "") { inputFile = countfile; }
        
        for(int i=0;i<numGroups;i++){  m->mothurRemove((inputFile + "." + toString(i) + ".temp")); }

        singleton = inputFile + ".extra.temp";
        ofstream remainingNames;
        m->openOutputFile(singleton, remainingNames);
        
        bool wroteExtra = false;
        
        ifstream bigNameFile;
        m->openInputFile(inputFile, bigNameFile);
        
        //grab header line 
        string headers = "";
        if (countfile != "") { headers = m->getline(bigNameFile); m->gobble(bigNameFile); }
        
        string name, nameList;
        while(!bigNameFile.eof()){
            bigNameFile >> name >> nameList;  
            m->getline(bigNameFile); m->gobble(bigNameFile); //extra getline is for rest of countfile line if groups are given.
            
            //did this sequence get assigned a group
            it = seqGroup.find(name);
            
            if (it != seqGroup.end()) {  
                m->openOutputFileAppend((inputFile + "." + toString(it->second) + ".temp"), outFile);
                outFile << name << '\t' << nameList << endl;
                outFile.close();
            }else{
                wroteExtra = true;
                remainingNames << name << '\t' << nameList << endl;
            }
        }
        bigNameFile.close();
        
		for(int i=0;i<numGroups;i++){
			string tempNameFile = inputFile + "." + toString(i) + ".temp";
			string tempDistFile = tempDistFiles[i];
            
            //if there are valid distances
            ifstream fileHandle;
            fileHandle.open(tempDistFile.c_str());
            if(fileHandle) 	{
                m->gobble(fileHandle);
                if (!fileHandle.eof()) {  //check
                    map<string, string> temp;
                    if (countfile != "") {
                        //add header
                        ofstream out;
                        string newtempNameFile = tempNameFile + "2";
                        m->openOutputFile(newtempNameFile, out);
                        out << "Representative_Sequence\ttotal" << endl;
                        out.close();
                        m->appendFiles(tempNameFile, newtempNameFile);
                        m->mothurRemove(tempNameFile);
                        m->renameFile(newtempNameFile, tempNameFile);
                    }
                    temp[tempDistFile] = tempNameFile;
                    dists.push_back(temp);
                }else{
                    ifstream in;
                    m->openInputFile(tempNameFile, in);
                    
                    while(!in.eof()) {
                        in >> name >> nameList;  m->gobble(in);
                        wroteExtra = true;
                        remainingNames << name << '\t' << nameList << endl;
                    }
                    in.close();
                    m->mothurRemove(tempNameFile);
                }
            }
            fileHandle.close();
		}
		
		remainingNames.close();
		
		if (!wroteExtra) { 
			m->mothurRemove(singleton);
			singleton = "none";
		}else if (countfile != "") {
            //add header
            ofstream out;
            string newtempNameFile = singleton + "2";
            m->openOutputFile(newtempNameFile, out);
            out << "Representative_Sequence\ttotal" << endl; 
            out.close();
            m->appendFiles(singleton, newtempNameFile);
            m->mothurRemove(singleton);
            m->renameFile(newtempNameFile, singleton);
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
		m->openInputFile(distFile, dFile);

		while(dFile){
			string seqA, seqB;
			float dist;

			dFile >> seqA >> seqB >> dist;
			
			if (m->control_pressed) {   dFile.close();  for(int i=0;i<numGroups;i++){	if(groups[i].size() > 0){  m->mothurRemove((distFile + "." + toString(i) + ".temp")); }  } return 0; }
					
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
			m->gobble(dFile);
		}
		dFile.close();
		
        vector<string> tempDistFiles;
		for (int i = 0; i < numGroups; i++) {
            string fileName = distFile + "." + toString(i) + ".temp";
            tempDistFiles.push_back(fileName);
			if (outputs[i] != "") {
				ofstream outFile;
				outFile.open(fileName.c_str(), ios::ate);
				outFile << outputs[i];
				outFile.close();
			}
		}
        
        map<string, int> seqGroup;
        for (int i = 0; i < groups.size(); i++) {
            for (set<string>::iterator itNames = groups[i].begin(); itNames != groups[i].end();) {
                seqGroup[*itNames] = i;
                groups[i].erase(itNames++);
            }
        }
        
		splitNames(seqGroup, numGroups, tempDistFiles);
				
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

