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
//column formatted distance file
SplitMatrix::SplitMatrix(string distfile, string name, string count, string tax, float c, string t){
	m = MothurOut::getInstance();
    
	distFile = distfile;
	namefile = name;
    taxFile = tax;
    countfile = count;
    
	method = t;
    cutoff = c;
    outputType = "distance";
}
/***********************************************************************/

SplitMatrix::SplitMatrix(string ffile, string name, string count, string tax, float c, float cu, string t, int p, bool cl, string output, string ot){
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
    outputType = ot;
}

/***********************************************************************/

int SplitMatrix::split(){
	try {
        
		if (method == "distance") {  
			splitDistance();
		}else if ((method == "classify") || (method == "fasta") || (method == "vsearch")) {
			splitClassify();
		}else {
			m->mothurOut("Unknown splitting method, aborting split.\n"); 
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
        util.readTax(taxFile, temp, true);
        
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
		util.openInputFile(fastafile, in);

		for (int i = 0; i < numGroups; i++) {  util.mothurRemove((fastafile + "." + toString(i) + ".temp")); }
	
		//parse fastafile
		while (!in.eof()) {
			Sequence query(in); util.gobble(in);
			if (query.getName() != "") {
		
				it = seqGroup.find(query.getName());
				
				//save names in case no namefile is given
				if ((namefile == "") && (countfile == "")) {  names.insert(query.getName()); }
			
				if (it != seqGroup.end()) { //not singleton
                    ofstream outFile;
                    util.openOutputFileAppend((fastafile + "." + toString(it->second) + ".temp"), outFile);
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
			m->mothurOut("ERROR: " + it->first + " is missing from your fastafile. This could happen if your taxonomy file is not unique and your fastafile is, or it could indicate an error.\n"); 
            error = true;
		}
		copyGroups.clear();
        
        if (error) { exit(1); }
        
        if (outputType == "distance") { //create distance matrices for each fasta file
            //process each distance file
            for (int i = 0; i < numGroups; i++) {
                
                string options = "";
                if (classic) { options += "fasta=" + (fastafile + "." + toString(i) + ".temp") + ", processors=" + toString(processors) + ", output=lt"; }
                else { options += "fasta=" + (fastafile + "." + toString(i) + ".temp") + ", processors=" + toString(processors) + ", cutoff=" + toString(distCutoff); }
                if (outputDir != "") { options += ", outputdir=" + outputDir; }
                
                m->mothurOut("/******************************************/\n");
                m->mothurOut("Running command: dist.seqs(" + options + ")\n");
                Command* command = new DistanceCommand(options);
                m->mothurOut("/******************************************/\n");
                
                command->execute();
                delete command;
                
                util.mothurRemove((fastafile + "." + toString(i) + ".temp"));
                
                //remove old names files just in case
                if (namefile != "") { util.mothurRemove((namefile + "." + toString(i) + ".temp")); }
                else { util.mothurRemove((countfile + "." + toString(i) + ".temp")); }
            }
        }
        
        vector<string> tempDistFiles;    
        for(int i=0;i<numGroups;i++){
            if (outputDir == "") { outputDir = util.hasPath(fastafile); }
            string tempDistFile = (fastafile + "." + toString(i) + ".temp");
            if (outputType == "distance") {
                if (classic) { tempDistFile =  outputDir + util.getRootName(util.getSimpleName((fastafile + "." + toString(i) + ".temp"))) + "phylip.dist";}
                else { tempDistFile = outputDir + util.getRootName(util.getSimpleName((fastafile + "." + toString(i) + ".temp"))) + "dist"; }
            }
            tempDistFiles.push_back(tempDistFile);
        }
        
        if (method == "vsearch")    {   splitNamesVsearch(seqGroup, numGroups, tempDistFiles);  }
        else                        {   splitNames(seqGroup, numGroups, tempDistFiles);     }
        
		if (m->getControl_pressed())	 {  for (int i = 0; i < dists.size(); i++) { util.mothurRemove((dists[i].begin()->first)); util.mothurRemove((dists[i].begin()->second)); } dists.clear(); }

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
		util.openInputFile(distFile, dFile);
		
		
		for (int i = 0; i < numGroups; i++) { //remove old temp files, just in case
			util.mothurRemove((distFile + "." + toString(i) + ".temp"));
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
			
			if (m->getControl_pressed()) { dFile.close(); for (int i = 0; i < numGroups; i++) { util.mothurRemove((distFile + "." + toString(i) + ".temp"));	} }
			
			dFile >> seqA >> seqB >> dist;  util.gobble(dFile);
			
			//if both sequences are in the same group then they are within the cutoff
			it = seqGroup.find(seqA);
			it2 = seqGroup.find(seqB);
			
			if ((it != seqGroup.end()) && (it2 != seqGroup.end())) { //they are both not singletons 
				if (it->second == it2->second) { //they are from the same group so add the distance
					if (numOutputs[it->second] > 30) {
						util.openOutputFileAppend((distFile + "." + toString(it->second) + ".temp"), outFile);
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
			util.mothurRemove((inputFile + "." + toString(i) + ".temp"));
			
			//write out any remaining buffers
			if (numOutputs[i] > 0) {
				util.openOutputFileAppend((distFile + "." + toString(i) + ".temp"), outFile);
				outFile << outputs[i];
				outFile.close();
				outputs[i] = "";
				numOutputs[i] = 0;
				validDistances[i] = true;
			}
		}
		
        splitNames(seqGroup, numGroups, tempDistFiles);
        
		if (m->getControl_pressed())	 {  
			for (int i = 0; i < dists.size(); i++) { 
				util.mothurRemove((dists[i].begin()->first));
				util.mothurRemove((dists[i].begin()->second));
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
//********************************************************************************************************************
int SplitMatrix::splitNames(map<string, int>& seqGroup, int numGroups, vector<string>& tempDistFiles){
	try {
        ofstream outFile;
        map<string, int>::iterator it;
        
        string inputFile = namefile;
        if (countfile != "") { inputFile = countfile; }
        
        for(int i=0;i<numGroups;i++){  util.mothurRemove((inputFile + "." + toString(i) + ".temp")); }

        singleton = inputFile + ".extra.temp";
        ofstream remainingNames;
        util.openOutputFile(singleton, remainingNames);
        
        bool wroteExtra = false;
        
        //grab header line
        string defaultCountTableHeaders = "";
        if (countfile != "") {
            map<string, int> ctMap;
            CountTable ct; ct.readTable(countfile, false, true); ctMap = ct.getNameMap();
            vector<string> headers = ct.getHardCodedHeaders();
            defaultCountTableHeaders = util.getStringFromVector(headers, "\t");
            for (map<string, int>::iterator itCtMap = ctMap.begin(); itCtMap != ctMap.end(); itCtMap++) {
                string name = itCtMap->first;
                int abundance = itCtMap->second;
                
                //did this sequence get assigned a group
                it = seqGroup.find(name);
                
                if (it != seqGroup.end()) {
                    util.openOutputFileAppend((inputFile + "." + toString(it->second) + ".temp"), outFile);
                    outFile << name << '\t' << abundance << endl;
                    outFile.close();
                }else{
                    wroteExtra = true;
                    remainingNames << name << '\t' << abundance << endl;
                    numSingleton++;
                }
            }
        }else {
            ifstream bigNameFile;
            util.openInputFile(inputFile, bigNameFile);
            
            string name, nameList;
            numSingleton = 0;
            while(!bigNameFile.eof()){
                bigNameFile >> name; util.gobble(bigNameFile); bigNameFile >> nameList; util.gobble(bigNameFile);
                
                //did this sequence get assigned a group
                it = seqGroup.find(name);
                
                if (it != seqGroup.end()) {
                    util.openOutputFileAppend((inputFile + "." + toString(it->second) + ".temp"), outFile);
                    outFile << name << '\t' << nameList << endl;
                    outFile.close();
                }else{
                    wroteExtra = true;
                    remainingNames << name << '\t' << nameList << endl;
                    numSingleton++;
                }
            }
            bigNameFile.close();
        }
        
        
		for(int i=0;i<numGroups;i++){
			string tempNameFile = inputFile + "." + toString(i) + ".temp";
			string tempDistFile = tempDistFiles[i];
            
            //if there are valid distances
            ifstream fileHandle;
            fileHandle.open(tempDistFile.c_str());
            
            bool removeDist = false;
            if(fileHandle) 	{
                util.gobble(fileHandle);
                if (!fileHandle.eof()) {  //check
                    map<string, string> temp;
                    if (countfile != "") {
                        //add header
                        ofstream out;
                        string newtempNameFile = tempNameFile + "2";
                        util.openOutputFile(newtempNameFile, out);
                        out << defaultCountTableHeaders << endl;
                        out.close();
                        util.appendFiles(tempNameFile, newtempNameFile);
                        util.mothurRemove(tempNameFile);
                        util.renameFile(newtempNameFile, tempNameFile);
                    }
                    temp[tempDistFile] = tempNameFile;
                    dists.push_back(temp);
                }else{
                    ifstream in;
                    util.openInputFile(tempNameFile, in);
                    
                    string name, nameList;
                    while(!in.eof()) {
                        in >> name >> nameList;  util.gobble(in);
                        wroteExtra = true;
                        remainingNames << name << '\t' << nameList << endl;
                        numSingleton++;
                    }
                    in.close();
                    util.mothurRemove(tempNameFile);
                    removeDist = true;
                }
            }
            fileHandle.close();
            if (removeDist) { util.mothurRemove(tempDistFile); }
		}
		
		remainingNames.close();
		
		if (!wroteExtra) { 
			util.mothurRemove(singleton);
			singleton = "none";
            numSingleton = 0;
		}else if (countfile != "") {
            //add header
            ofstream out;
            string newtempNameFile = singleton + "2";
            util.openOutputFile(newtempNameFile, out);
            out << defaultCountTableHeaders << endl;
            out.close();
            util.appendFiles(singleton, newtempNameFile);
            util.mothurRemove(singleton);
            util.renameFile(newtempNameFile, singleton);
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitNames");
		exit(1);
	}
}
//********************************************************************************************************************
int SplitMatrix::splitNamesVsearch(map<string, int>& seqGroup, int numGroups, vector<string>& tempDistFiles){
    try {
        ofstream outFile;
        map<string, int>::iterator it;
        
        string inputFile = namefile;
        if (countfile != "") { inputFile = countfile; }
        
        for(int i=0;i<numGroups;i++){  util.mothurRemove((inputFile + "." + toString(i) + ".temp")); }
        
        singleton = inputFile + ".extra.temp";
        ofstream remainingNames;
        util.openOutputFile(singleton, remainingNames);
        
        bool wroteExtra = false;
        string errorMessage = "name";
        
        //grab header line
        string name, nameList;
        string defaultCountTableHeaders = "";
        numSingleton = 0;
        string headers = ""; bool hasGroups = false;
        if (countfile != "") {
            errorMessage = "count"; map<string, int> ctMap;
            CountTable ct; ct.readTable(countfile, true, true); ctMap = ct.getNameMap();
            hasGroups = ct.hasGroupInfo();
            vector<string> headers = ct.getHardCodedHeaders();
            defaultCountTableHeaders = util.getStringFromVector(headers, "\t");
            
            for (map<string, int>::iterator itCtMap = ctMap.begin(); itCtMap != ctMap.end(); itCtMap++) {
                string name = itCtMap->first;
                int abundance = itCtMap->second;
                
                //did this sequence get assigned a group
                it = seqGroup.find(name);
                
                if (it != seqGroup.end()) {
                    util.openOutputFileAppend((inputFile + "." + toString(it->second) + ".temp"), outFile);
                    outFile << name << '\t' << abundance << endl;
                    outFile.close();
                }else{
                    wroteExtra = true;
                    remainingNames << name << '\t' << abundance << endl;
                    numSingleton++;
                }
                
            }
        }else {
            ifstream bigNameFile;
            util.openInputFile(inputFile, bigNameFile);
            
            while(!bigNameFile.eof()){
                bigNameFile >> name >> nameList; util.gobble(bigNameFile);
                if (hasGroups) { util.getline(bigNameFile); util.gobble(bigNameFile);  }
                
                //did this sequence get assigned a group
                it = seqGroup.find(name);
                
                if (it != seqGroup.end()) {
                    util.openOutputFileAppend((inputFile + "." + toString(it->second) + ".temp"), outFile);
                    outFile << name << '\t' << nameList << endl;
                    outFile.close();
                }else{
                    wroteExtra = true;
                    remainingNames << name << '\t' << nameList << endl;
                    numSingleton++;
                }
            }
            bigNameFile.close();

        }

        for(int i=0;i<numGroups;i++){
            string tempNameFile = inputFile + "." + toString(i) + ".temp";
            string tempDistFile = tempDistFiles[i];
            
            //if there are valid distances
            ifstream fileHandle;
            fileHandle.open(tempDistFile.c_str());
            if(fileHandle) 	{
                util.gobble(fileHandle);
                if (!fileHandle.eof()) {  //check
                    map<string, string> temp;
                    if (countfile != "") {
                        //add header
                        ofstream out;
                        string newtempNameFile = tempNameFile + "2";
                        util.openOutputFile(newtempNameFile, out);
                        out << defaultCountTableHeaders << endl;
                        out.close();
                        util.appendFiles(tempNameFile, newtempNameFile);
                        util.mothurRemove(tempNameFile);
                        util.renameFile(newtempNameFile, tempNameFile);
                    }
                    temp[tempDistFile] = tempNameFile;
                    dists.push_back(temp);
                }else{
                    ifstream in;
                    util.openInputFile(tempNameFile, in);
                    
                    while(!in.eof()) {
                        in >> name >> nameList;  util.gobble(in);
                        wroteExtra = true;
                        remainingNames << name << '\t' << nameList << endl;
                        numSingleton++;
                    }
                    in.close();
                    util.mothurRemove(tempNameFile);
                }
            }
            fileHandle.close();
        }
        
        remainingNames.close();
        
        if (!wroteExtra) {
            util.mothurRemove(singleton);
            singleton = "none";
            numSingleton = 0;
        }else if (countfile != "") {
            //add header
            ofstream out;
            string newtempNameFile = singleton + "2";
            util.openOutputFile(newtempNameFile, out);
            out << defaultCountTableHeaders << endl;
            out.close();
            util.appendFiles(singleton, newtempNameFile);
            util.mothurRemove(singleton);
            util.renameFile(newtempNameFile, singleton);
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitMatrix", "splitNamesTax");
        exit(1);
    }
}
//********************************************************************************************************************
int SplitMatrix::splitDistance(){
	try {
        Utils util;
        vector<string> nameMap;
        map<string, long long> nameAssignment;
        if (namefile != "") {
            map<string, string> namefileMap;
            util.readNames(namefile, namefileMap);
            for (map<string, string>::iterator it = namefileMap.begin(); it != namefileMap.end(); it++) {
                nameAssignment[it->first] = 1; //value of 1 reset below, just a placeholder
            }
        }
        else  {
            CountTable ct; ct.readTable(countfile, false, true);
            map<string, int> temp = ct.getNameMap();
            for (map<string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
        }
        long long count = 0;
        for (map<string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
            it->second = count; count++;
            nameMap.push_back(it->first);
        }
        
        int numSeqs = nameMap.size();
        
        //seqsseqsGroupAssignment[0] = group assignment for seq 0
        vector<long long> seqsGroupAssignment(numSeqs, -1); //assign all seqs no group
        
        //when we merge groups, rather than reassigning all the seqs in that group, let's reassign the group
        //mergedGroups[0] = group assignment for group 0.
        vector<int> mergedGroups;
        
        vector<string> groupFiles; //filenames for each group
        vector<string> groupOutputs; //buffered file outputs for each group
        vector<int> groupOutputCounts; //count to empty buffer for each group
		
		int numGroups = 0;
    
		ifstream dFile;
		util.openInputFile(distFile, dFile);

		while(dFile){
			string seqA, seqB;
			float dist;

            dFile >> seqA >> seqB >> dist; util.gobble(dFile);
			
			if (m->getControl_pressed()) {  break; }
					
			if(dist <= cutoff){
                
                map<string,long long>::iterator itA = nameAssignment.find(seqA);
                map<string,long long>::iterator itB = nameAssignment.find(seqB);
                
                if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + seqA + "' was not found in the name or count file, please correct\n"); exit(1);  }
                if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + seqB + "' was not found in the name or count file, please correct\n"); exit(1);  }

                long long indexA = (itA->second);
                long long indexB = (itB->second);
                
                int groupIDA = findRootGroup(mergedGroups, seqsGroupAssignment[indexA]);
                int groupIDB = findRootGroup(mergedGroups, seqsGroupAssignment[indexB]);
                int groupID = -1;
                
                if(groupIDA != -1 && groupIDB != -1){ //both are already assigned to a group, so merge the groups. set to groupIDA
                    
                    int thisGroup = mergedGroups[groupIDA];
                    int thatGroup = mergedGroups[groupIDB];
                    
                    //merge files and save in groupIDA
                    util.appendFiles(groupFiles[thatGroup], groupFiles[thisGroup]);
                    ofstream out; util.openOutputFile(groupFiles[thatGroup], out); out.close(); //clear file
                    
                    //append outputs
                    groupOutputs[thisGroup] += groupOutputs[thatGroup];
                    groupOutputs[thatGroup] = "";
                    
                    //update output counts
                    groupOutputCounts[thisGroup] += groupOutputCounts[thatGroup];
                    groupOutputCounts[thatGroup] = 0;
                    
                    //set group id
                    groupID = thisGroup;
                    mergedGroups[groupIDB] = thisGroup;
            
                }else if(groupIDA != -1 && groupIDB == -1){ //seqA has a group and seqB is unassigned
                    groupID = mergedGroups[groupIDA]; //assign seqB to seqA's group
                    
                }else if(groupIDA == -1 && groupIDB != -1) {  //seqB has a group and seqA is unassigned
                    
                    groupID = mergedGroups[groupIDB]; //assign seqA to seqB's group
                    
                }else {  //both seqs have no group
                    groupID = numGroups; //we need a new group
                    mergedGroups.push_back(numGroups); //assign group to merge with self
                    string fileName = distFile + "." + toString(numGroups) + ".temp";
                    ofstream out; util.openOutputFile(fileName, out); out.close(); //clear file
                    groupFiles.push_back(fileName);
                    groupOutputs.push_back("");
                    groupOutputCounts.push_back(0);
                    numGroups++;
                }
                
                string output = seqA + '\t' + seqB +'\t' + toString(dist) + '\n';
                
                seqsGroupAssignment[indexA] = groupID;
                seqsGroupAssignment[indexB] = groupID;
                groupOutputs[groupID] += output;
                groupOutputCounts[groupID]++;
                
                if (groupOutputCounts[groupID] > 10000) {
                    ofstream out; util.openOutputFileAppend(groupFiles[groupID], out);
                    out << groupOutputs[groupID]; out.close();
                    groupOutputCounts[groupID] = 0;
                    groupOutputs[groupID] = "";
                }
			}
			
		}
		dFile.close();
        
        vector<string> tempDistFiles;
        for (int i = 0; i < numGroups; i++) {
            if (groupOutputCounts[i] != 0) { //write remaning buffer for group
                ofstream out; util.openOutputFileAppend(groupFiles[i], out);
                out << groupOutputs[i]; out.close();
                groupOutputCounts[i] = 0;
                groupOutputs[i] = "";
                tempDistFiles.push_back(groupFiles[i]);
            }else {
                //check for blank group - happens when groups merge
                if (util.isBlank(groupFiles[i])) {
                    util.mothurRemove(groupFiles[i]);
                }else {
                    tempDistFiles.push_back(groupFiles[i]);
                }
            }
        }
        
        map<string, int> seqGroup;
        for (long long i = 0; i < seqsGroupAssignment.size(); i++) {
            
            if (m->getControl_pressed()) {  break;  }
            
            if (seqsGroupAssignment[i] != -1) { //you have a distance below the cutoff
                int group = findRootGroup(mergedGroups, seqsGroupAssignment[i]);
                
                seqGroup[nameMap[i]] = group;
            }
        }
    
		splitNames(seqGroup, numGroups, tempDistFiles);
				
		return 0;			
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitDistance");
		exit(1);
	}
}
/***********************************************************************/
int SplitMatrix::findRootGroup(vector<int>& mergedGroups, int pos){
    try {
        
        //if unassigned
        if (pos == -1) { return pos; }
        
        //if the mergedGroups[10] = 10 then you are at the root
        //if mergedGroups[10] != 10, find parent group
        //mergedGroups[10] = 5, then check mergeGroups[5].
        
        int rootGroup = mergedGroups[pos];
        
        if (rootGroup != pos) { //you need to look at your parent
            rootGroup = findRootGroup(mergedGroups, rootGroup);
            mergedGroups[pos] = rootGroup;
        }//else you are the root
        
        return rootGroup;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitMatrix", "findRootGroup");
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

