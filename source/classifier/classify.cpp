/*
 *  classify.cpp
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "classify.h"
#include "sequence.hpp"
#include "kmerdb.hpp"
#include "suffixdb.hpp"
#include "distancedb.hpp"
#include "optidb.hpp"


/**************************************************************************************************/
void Classify::generateDatabaseAndNames(string tfile, string tempFile, string method, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch, string version)  {
	try {
        Utils util;
        
        int numSeqs = 0; maxLevel = 0;
		taxFile = tfile;
        templateFile = tempFile;

        long start = time(nullptr); m->mothurOut("Generating search database...    "); cout.flush();
        
        //need to know number of template seqs for suffixdb
        if (method == "suffix") {
            ifstream inFASTA; util.openInputFile(tempFile, inFASTA);
            util.getNumSeqs(inFASTA, numSeqs);
            inFASTA.close();
        }
        
        bool needToGenerate = true;
        string dBName;
        if(method == "kmer")			{
            database = new KmerDB(tempFile, kmerSize);
            
            dBName = tempFile.substr(0,tempFile.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
            ifstream kmerFileTest(dBName.c_str());
            if(kmerFileTest){
                string line = util.getline(kmerFileTest);
                bool GoodFile = util.checkReleaseVersion(line, version); kmerFileTest.close();
                int shortcutTimeStamp = util.getTimeStamp(dBName);
                int referenceTimeStamp = util.getTimeStamp(tempFile);
                
                //if the shortcut file is older then the reference file, remake shortcut file
                if (shortcutTimeStamp < referenceTimeStamp) {  GoodFile = false;  }

                if (GoodFile) {  needToGenerate = false;	}
            }
        }
        else if(method == "suffix")		{	database = new SuffixDB(numSeqs);								}
        else if(method == "distance")	{	database = new DistanceDB();	}
        else {
            m->mothurOut(method + " is not a valid search option. I will run the command using kmer, ksize=8.\n");
            database = new KmerDB(tempFile, 8);
        }
        
        if (!m->getControl_pressed()) {
            if (needToGenerate) {
                ifstream fastaFile; util.openInputFile(tempFile, fastaFile);
                
                while (!fastaFile.eof()) {
                    Sequence temp(fastaFile); util.gobble(fastaFile);
                    
                    names.push_back(temp.getName());
                    
                    database->addSequence(temp);
                }
                fastaFile.close();
                
                if ((method == "kmer") && (!shortcuts)) {;} //don't print
                else {database->generateDB(); }
                
            }else if ((method == "kmer") && (!needToGenerate)) {
                ifstream FileTest(dBName.c_str());
                database->readDB(FileTest);
                
                ifstream fastaFile; util.openInputFile(tempFile, fastaFile);
                
                while (!fastaFile.eof()) {
                    Sequence temp(fastaFile); util.gobble(fastaFile);
                    
                    names.push_back(temp.getName());
                }
                fastaFile.close();
            }
            database->setNumSeqs(names.size());
            
            m->mothurOut("DONE.\nIt took " + toString(time(nullptr) - start) + " seconds generate search database.\n");
            
            readTaxonomy(taxFile);
            
            //sanity check
            bool okay = phyloTree->ErrorCheck(names);
            
            if (!okay) { m->setControl_pressed(true); }
        }
	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "generateDatabaseAndNames");
		exit(1);
	}
}
/**************************************************************************************************/
Classify::Classify() {		m = MothurOut::getInstance();	database = nullptr;	phyloTree=nullptr;   maxLevel = 0; }
/**************************************************************************************************/

int Classify::readTaxonomy(string file) {
	try {
		phyloTree = new PhyloTree();
		string name, taxInfo;
		
		m->mothurOut("\nReading in the " + file + " taxonomy...\t");	cout.flush();
        if (m->getDebug()) { m->mothurOut("[DEBUG]: Taxonomies read in...\n"); }
        
        taxonomy.clear();
        
        Utils util; util.readTax(file, taxonomy, true);
        
        for (map<string, string>::iterator itTax = taxonomy.begin(); itTax != taxonomy.end(); itTax++) {
            phyloTree->addSeqToTree(itTax->first, itTax->second);
            if (m->getControl_pressed()) { break; }
        }
       
		phyloTree->assignHeirarchyIDs(0);
        
        maxLevel = phyloTree->getMaxLevel();
		
		phyloTree->setUp(file);

		m->mothurOut("DONE.\n"); cout.flush();
		
		return phyloTree->getNumSeqs();
	
	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "readTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

vector<string> Classify::parseTax(string tax) {
	try {
		vector<string> taxons;
		Utils util; util.splitAtChar(tax, taxons, ';');
        return taxons;
	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "parseTax");
		exit(1);
	}
}
/**************************************************************************************************/

double Classify::getLogExpSum(vector<double> probabilities, int& maxIndex){
	try {
        //	http://jblevins.org/notes/log-sum-exp
        
        double maxProb = probabilities[0];
        maxIndex = 0;
        
        int numProbs = (int)probabilities.size();
        
        for(int i=1;i<numProbs;i++){
            if(probabilities[i] >= maxProb){
                maxProb = probabilities[i];
                maxIndex = i;
            }
        }
        
        double probSum = 0.0000;
        
        for(int i=0;i<numProbs;i++){
            probSum += exp(probabilities[i] - maxProb);		
        }
        
        probSum = log(probSum) + maxProb;
        
        return probSum;
	}
	catch(exception& e) {
		m->errorOut(e, "Classify", "getLogExpSum");
		exit(1);
	}
}
/**************************************************************************************************/
bool Classify::checkReleaseDate(vector<ifstream*>& files, string version) {
    try {
        Utils util;
        bool good = true;
        
        vector<string> versionVector;
        util.splitAtChar(version, versionVector, '.');
        
        for (int i = 0; i < files.size(); i++) {
            string line = util.getline(*files[i]);
            
            if (line[0] != '#') { good = false; break; } //shortcut files from before we added this check
            else {
                line = line.substr(1);
                
                vector<string> linesVector; util.splitAtChar(line, linesVector, '.');
                
                if (versionVector.size() != linesVector.size()) { good = false; break; }
                else {
                    for (int j = 0; j < versionVector.size(); j++) {
                        int num1, num2;
                        convert(versionVector[j], num1);
                        convert(linesVector[j], num2);
                        
                        //if mothurs version is newer than this files version, then we want to remake it
                        if (num1 > num2) {  good = false; break;  }
                    }
                }
                if (!good) { break; }
            }
        }
        
        if (!good) {  for (int i = 0; i < files.size(); i++) { files[i]->close(); }  }
        else { for (int i = 0; i < files.size(); i++) { files[i]->seekg(0); }  }
        
        return good;
    }
    catch(exception& e) {
        m->errorOut(e, "Classify", "checkReleaseDate");
        exit(1);
    }
}
/**************************************************************************************************/

