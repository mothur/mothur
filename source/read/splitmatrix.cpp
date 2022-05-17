/*
 *  splitmatrix.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 */


#include "splitmatrix.h"
#include "phylotree.h"
#include "distancecommand.h"
#include "pairwiseseqscommand.h"
#include "seqsummarycommand.h"
#include "getseqscommand.h"
#include "removeseqscommand.h"

 /***********************************************************************/
SplitMatrix::SplitMatrix(string ffile, string name, string count, string tax, float c, float cu, int p, bool cl, string output, bool v){
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
    usingVsearchToCLuster = v;
    
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
		
        if (usingVsearchToCLuster)  { createFastaFilesFromTax(seqGroups, taxGroupNames);        }
        else                        {  createDistanceFilesFromTax(seqGroups, taxGroupNames);    }
		
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
        ofstream outNonSingleton; util.openOutputFile(nonSingletonsFile, outNonSingleton);
        
        ifstream inFASTA; util.openInputFile(fastafile, inFASTA);
        SequenceDB fullDB(inFASTA); inFASTA.close();
        
        if (m->getDebug()) { for (int i = 0; i < numGroups; i++) { m->mothurOut("[DEBUG]: Number of unique sequences for group " + groupNames[i] + " (" + toString(i+1) + " of " + toString(numGroups) + "): " + toString(seqGroups[i].size()) + "\n\n"); } }
        
        //process each group
        for (int i = 0; i < numGroups; i++) {
            
            if (m->getControl_pressed()) { outNonSingleton.close(); util.mothurRemove(nonSingletonsFile); for (int i = 0; i < dists.size(); i++) { util.mothurRemove((dists[i].begin()->first)); util.mothurRemove((dists[i].begin()->second)); } dists.clear(); return 0; }
            
            unordered_set<string> thisGroupsNames = util.mothurConvert(seqGroups[i]);
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Selecting sequences for group " + groupNames[i] + " (" + toString(i+1) + " of " + toString(numGroups) + ")\nNumber of unique sequences: " + toString(seqGroups[i].size()) + "\n\n");
            
            string outName = "";
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
                outName = thisOutputDir + util.getRootName(util.getSimpleName(namefile)) + toString(i) + ".name.temp";
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
                outName = thisOutputDir + util.getRootName(util.getSimpleName(countfile)) + toString(i) + ".count.temp";
            }
            
            pair<string,string> dupsFile(namefile, outName); string dupsFormat = "name";
            if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
            
            Command* getCommand = new GetSeqsCommand(thisGroupsNames, nullStringPair, nullStringPair, dupsFile, dupsFormat);
            
            delete getCommand;
            
            StorageDatabase* thisDB;
            
            m->mothurOut("\nCalculating distances for group " + groupNames[i] + " (" + toString(i+1) + " of " + toString(numGroups) + "):\n");
            thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
            string outputFileRoot = thisOutputDir + util.getRootName(util.getSimpleName(fastafile)) + toString(i) + ".";
            
            string outputformat = "column"; if (classic) { outputformat = "lt"; }
            
            Command* command;
            vector< vector< int > > kmerDB; vector< int > lengths;
            
            if (fullDB.sameLength()) {
                thisDB = new SequenceDB(fullDB, thisGroupsNames);
                command = new DistanceCommand(thisDB, outputFileRoot, distCutoff, outputformat, processors);
            }
            else {
                thisDB = new SequenceDB(fullDB, thisGroupsNames, 7, kmerDB, lengths);
                command = new PairwiseSeqsCommand(thisDB, kmerDB, lengths, outputFileRoot, distCutoff, outputformat, processors);
            }
            
            map<string, vector<string> > filenames = command->getOutputFiles();
            
            string thisDistanceFile = "";
            if (classic) { thisDistanceFile = filenames["phylip"][0]; }
            else { thisDistanceFile = filenames["column"][0]; }
            
            delete command;
            
            m->mothurOut("/******************************************/\n");
            
            if (!util.isBlank(thisDistanceFile)) {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
                string outDist = thisOutputDir + util.getRootName(util.getSimpleName(fastafile)) + toString(i) + ".disttemp";
                util.renameFile(thisDistanceFile, outDist);
                
                map<string, string> thisFilePair; thisFilePair[outDist] = outName;
                dists.push_back(thisFilePair);
                
                for (int j = 0; j < thisDB->getNumSeqs(); j++) { outNonSingleton << thisDB->getSeq(j).getName() << endl; }
            }else {
                util.mothurRemove(thisDistanceFile);
                util.mothurRemove(outName);
            }
            
            delete thisDB;
        }
        outNonSingleton.close();
        
        
        if (!util.isBlank(nonSingletonsFile)) { //there are non singletons, so remove them to find the singletons
            //get singletons
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
                singleton = thisOutputDir + util.getRootName(util.getSimpleName(namefile)) + "singletons.temp";
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
                singleton = thisOutputDir + util.getRootName(util.getSimpleName(countfile)) + "singletons.temp";
            }
            
            pair<string,string> dupsFile(namefile, singleton); string dupsFormat = "name";
            if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Finding singletons (ignore 'Removing group' messages):\n\nRunning command: remove.seqs()\n");
            
            Command* removeCommand = new RemoveSeqsCommand(nonSingletonsFile, dupsFile, dupsFormat);
            
            delete removeCommand;
            
            m->mothurOut("/******************************************/\n");
            
        }else { //every seqs is a singleton
            if (namefile != "") { singleton = namefile; }
            else if (countfile != "") { singleton = countfile; }
        }
        
        if (util.isBlank(singleton)) {
            util.mothurRemove(singleton);
            singleton = "none";
        }
        
        util.mothurRemove(nonSingletonsFile);
        
        if (m->getControl_pressed())     {  for (int i = 0; i < dists.size(); i++) { util.mothurRemove((dists[i].begin()->first)); util.mothurRemove((dists[i].begin()->second)); } dists.clear(); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
        exit(1);
    }
}
/***********************************************************************/
int SplitMatrix::createFastaFilesFromTax(vector<vector<string> >& seqGroups, vector<string> groupNames){
    try {
        
        int numGroups = seqGroups.size();
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        singleton = "none";
          
        ifstream inFASTA; util.openInputFile(fastafile, inFASTA);
        SequenceDB fullDB(inFASTA); inFASTA.close();
        
        if (!fullDB.sameLength()) {
            m->mothurOut("[ERROR]: Cannot cluster using vsearch with unaligned sequences, please correct.\n\n"); m->setControl_pressed(true);
        }
        
        if (m->getDebug()) { for (int i = 0; i < numGroups; i++) { m->mothurOut("[DEBUG]: Number of unique sequences for group " + groupNames[i] + " (" + toString(i+1) + " of " + toString(numGroups) + "): " + toString(seqGroups[i].size()) + "\n\n"); } }
        
        //process each group
        for (int i = 0; i < numGroups; i++) {
            
            if (m->getControl_pressed()) {  for (int i = 0; i < dists.size(); i++) { util.mothurRemove((dists[i].begin()->first)); util.mothurRemove((dists[i].begin()->second)); } dists.clear(); return 0; }
            
            unordered_set<string> thisGroupsNames = util.mothurConvert(seqGroups[i]);
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Selecting sequences for group " + groupNames[i] + " (" + toString(i+1) + " of " + toString(numGroups) + ")\nNumber of unique sequences: " + toString(seqGroups[i].size()) + "\n\n");
            
            string outName = "";
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
                outName = thisOutputDir + util.getRootName(util.getSimpleName(namefile)) + toString(i) + ".name.temp";
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
                outName = thisOutputDir + util.getRootName(util.getSimpleName(countfile)) + toString(i) + ".count.temp";
            }
            
            pair<string, string> dupsFile(namefile, outName); string dupsFormat = "name";
            if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
            
            Command* getCommand = new GetSeqsCommand(thisGroupsNames, nullStringPair, nullStringPair, dupsFile, dupsFormat);
            
            delete getCommand;
            
            StorageDatabase* thisDB;
            
            thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
            string outFasta = thisOutputDir + util.getRootName(util.getSimpleName(fastafile)) + toString(i) + ".fastatemp";
            thisDB = new SequenceDB(fullDB, thisGroupsNames);
            thisDB->print(outFasta);
            delete thisDB;
            
            map<string, string> thisFilePair; thisFilePair[outFasta] = outName;
            dists.push_back(thisFilePair);
        }
        
        if (m->getControl_pressed())     {  for (int i = 0; i < dists.size(); i++) { util.mothurRemove((dists[i].begin()->first)); util.mothurRemove((dists[i].begin()->second)); } dists.clear(); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
        exit(1);
    }
}
/********************************************************************************************************************/
//sorts biggest to smallest
inline bool compareFileSizes(map<string, string> left, map<string, string> right){
	
	FILE * pFile;
	long leftsize = 0;
		
	//get num bytes in file
	string filename = left.begin()->first;
	pFile = fopen (filename.c_str(),"rb");
	string error = "Error opening " + filename;
	if (pFile==nullptr) perror (error.c_str());
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
	if (pFile2==nullptr) perror (error.c_str());
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
//********************************************************************************************************************/


