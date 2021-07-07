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
#include "getseqscommand.h"
#include "removeseqscommand.h"

/***********************************************************************/

SplitMatrix::SplitMatrix(string ffile, string name, string count, string tax, float c, float cu, int p, bool cl, string output){
	m = MothurOut::getInstance();
	fastafile = ffile;
	namefile = name;
    countfile = count;
	taxFile = tax;
	cutoff = c;  //tax level cutoff
	distCutoff = cu; //for fasta method if you are creating distance matrix you need a cutoff for that
	processors = p;
    classic = cl;
	outputDir = output;
    
    splitClassify();
   
}
/***********************************************************************/
void SplitMatrix::splitClassify(){
	try {
		cutoff = int(cutoff);
		
        map<string, string> temp; util.readTax(taxFile, temp, true);
        
        PhyloTree phylo;
        for (map<string, string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {
            
            if (m->getControl_pressed()) { return; }
            
            phylo.addSeqToTree(itTemp->first, itTemp->second);
            temp.erase(itTemp++);
        }
		
		phylo.assignHeirarchyIDs(0);

		//make sure the cutoff is not greater than maxlevel
		if (cutoff > phylo.getMaxLevel()) { m->mothurOut("splitcutoff is greater than the longest taxonomy, using " + toString(phylo.getMaxLevel())); m->mothurOutEndLine(); cutoff = phylo.getMaxLevel(); }
	
        vector<vector<string> > seqGroups; //seqFroups[0] -> vector of string containing names of seqs assigned to group 0
        vector<string> taxGroupNames;
		//for each node in tree
		for (int i = 0; i < phylo.getNumNodes(); i++) {
		
            if (m->getControl_pressed()) { return; }
            
			//is this node within the cutoff
			TaxNode taxon = phylo.get(i);
	
			if (taxon.level == cutoff) {//if yes, then create group containing this nodes sequences
				if (taxon.accessions.size() > 1) { //if this taxon just has one seq its a singleton
                    vector<string> thisGroupsSeqs;
                    for (int j = 0; j < taxon.accessions.size(); j++) {
                        thisGroupsSeqs.push_back(taxon.accessions[j]);
					}
                    seqGroups.push_back(thisGroupsSeqs);
                    taxGroupNames.push_back(taxon.name);
				}
			}
		}
		
        createDistanceFilesFromTax(seqGroups, taxGroupNames);
		
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "splitClassify");
		exit(1);
	}
}
/***********************************************************************/
int SplitMatrix::createDistanceFilesFromTax(vector<vector<string> >& seqGroups, vector<string> groupNames){
	try {
        
        int numGroups = seqGroups.size();
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        string nonSingletonsFile = thisOutputDir + util.getRootName(util.getSimpleName(fastafile)) + "nonsingleton.accnos";
        
        thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        string tempAccnos = thisOutputDir + util.getRootName(util.getSimpleName(fastafile)) + "accnos.temp";
        
        //process each group
        for (int i = 0; i < numGroups; i++) {
            
            util.printAccnos(tempAccnos, seqGroups[i]);
            
            //run get.seqs to select fasta and name or count files
            string options = "fasta=" + fastafile + ", accnos=" + tempAccnos;
            if (countfile != "") { options += ", count=" + countfile; }
            else if (namefile != "") { options += ", name=" + namefile, + ", dups=T"; }
            
            if (outputDir != "") { options += ", outputdir=" + outputDir; }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Selecting sequences for group " + groupNames[i] + ":\n\nRunning command: get.seqs(" + options + ")\n");
            
            Command* getCommand = new GetSeqsCommand(options);
            getCommand->execute();
            
            map<string, vector<string> > filenames = getCommand->getOutputFiles();
            
            string thisGroupsFastaFile = filenames["fasta"][0];
            string outName = "";
            
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
                outName = thisOutputDir + util.getRootName(util.getSimpleName(namefile)) + toString(i) + ".name.temp";
                util.renameFile(filenames["name"][0], outName);
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
                outName = thisOutputDir + util.getRootName(util.getSimpleName(countfile)) + toString(i) + ".count.temp";
                util.renameFile(filenames["count"][0], outName);
            }
            delete getCommand;
            
            m->mothurOut("/******************************************/\n");
            
            //run dist.seqs to calc the distances for each group
            if (classic) { options = "fasta=" + thisGroupsFastaFile + ", processors=" + toString(processors) + ", output=lt"; }
            else { options = "fasta=" + thisGroupsFastaFile + ", processors=" + toString(processors) + ", cutoff=" + toString(distCutoff); }
            if (outputDir != "") { options += ", outputdir=" + outputDir; }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Calculating distances for group " + groupNames[i] + ":\n\nRunning command: dist.seqs(" + options + ")\n");
            Command* command = new DistanceCommand(options);
            
            m->mothurOut("/******************************************/\n");
            
            command->execute();
            
            filenames = command->getOutputFiles();
            
            string thisDistanceFile = "";
            if (classic) { thisDistanceFile = filenames["phylip"][0]; }
            else { thisDistanceFile = filenames["column"][0]; }
            
            delete command;
            
            util.mothurRemove(thisGroupsFastaFile);
            
            if (!util.isBlank(thisDistanceFile)) {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
                string outDist = thisOutputDir + util.getRootName(util.getSimpleName(fastafile)) + toString(i) + ".disttemp";
                util.renameFile(thisDistanceFile, outDist);
                
                map<string, string> thisFilePair; thisFilePair[outDist] = outName;
                dists.push_back(thisFilePair);
                
                util.appendFiles(tempAccnos, nonSingletonsFile);
            }else {
                util.mothurRemove(thisDistanceFile);
                util.mothurRemove(outName);
            }
        }
       
        if (!util.isBlank(nonSingletonsFile)) { //there are non singletons, so remove them to find the singletons
            //get singletons
            //run remove.seqs to remove nonsingletons name or count files
            string options = "accnos=" + nonSingletonsFile;
            if (countfile != "") { options += ", count=" + countfile; }
            else if (namefile != "") { options += ", name=" + namefile, + ", dups=T"; }
            
            if (outputDir != "") { options += ", outputdir=" + outputDir; }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Finding singletons (ignore 'Removing group' messages):\n\nRunning command: remove.seqs(" + options + ")\n");
            
            Command* removeCommand = new RemoveSeqsCommand(options);
            removeCommand->execute();
            
            map<string, vector<string> > filenames = removeCommand->getOutputFiles();
            
            string outName = "";
            
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
                singleton = thisOutputDir + util.getRootName(util.getSimpleName(namefile)) + "singletons.temp";
                util.renameFile(filenames["name"][0], singleton);
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
                singleton = thisOutputDir + util.getRootName(util.getSimpleName(countfile)) + "singletons.temp";
                util.renameFile(filenames["count"][0], singleton);
            }
            delete removeCommand;
            
            m->mothurOut("/******************************************/\n");
            
        }else { //every seqs is a singleton
            if (namefile != "") { singleton = namefile; }
            else if (countfile != "") { singleton = countfile; }
        }
        
        util.mothurRemove(nonSingletonsFile);
        
        if (util.isBlank(singleton)) {
            util.mothurRemove(singleton);
            singleton = "none";
        }
        
		if (m->getControl_pressed())	 {  for (int i = 0; i < dists.size(); i++) { util.mothurRemove((dists[i].begin()->first)); util.mothurRemove((dists[i].begin()->second)); } dists.clear(); }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
		exit(1);
	}
}
/********************************************************************************************************************
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
}*/
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

