/*
 *  bayesian.cpp
 *  Mothur
 *
 *  Created by westcott on 11/3/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "bayesian.h"
#include "kmer.hpp"
#include "phylosummary.h"

/**************************************************************************************************/
Bayesian::Bayesian(string txfile, string tempFile, string method, int ksize, int cutoff, int i, int tid, bool f, bool sh, string version) :
Classify(), kmerSize(ksize), confidenceThreshold(cutoff), iters(i) {
	try {
		
		threadID = tid;
		flip = f;
        shortcuts = sh;
        string baseName = tempFile;
        string baseTName = txfile;
        Utils util;
        
        /************calculate the probablity that each word will be in a specific taxonomy*************/
        string tfileroot = util.getFullPathName(baseTName.substr(0,baseTName.find_last_of(".")+1));
        string tempfileroot = util.getRootName(util.getSimpleName(baseName));
        string phyloTreeName = tfileroot + "tree.train";
        string phyloTreeSumName = tfileroot + "tree.sum";
        string probFileName = tfileroot + tempfileroot + char('0'+ kmerSize) + "mer.prob";
        string probFileName2 = tfileroot + tempfileroot + char('0'+ kmerSize) + "mer.numNonZero";
        
        ofstream out;
        ofstream out2;
        
        vector<ifstream*> files;
        ifstream* phyloTreeTest = new ifstream(phyloTreeName.c_str()); files.push_back(phyloTreeTest);
		ifstream* probFileTest2 = new ifstream(probFileName2.c_str()); files.push_back(probFileTest2);
		ifstream* probFileTest = new ifstream(probFileName.c_str());   files.push_back(probFileTest);
		ifstream* probFileTest3 = new ifstream(phyloTreeSumName.c_str()); files.push_back(probFileTest3);
		
		long start = time(nullptr);
		
		//if they are there make sure they were created after this release date
		bool FilesGood = false;
		if(probFileTest && probFileTest2 && phyloTreeTest && probFileTest3){ FilesGood = checkReleaseDate(files, version); }

		if(probFileTest && probFileTest2 && phyloTreeTest && probFileTest3 && FilesGood){
			
			m->mothurOut("Reading template taxonomy...     "); cout.flush();
			
			phyloTree = new PhyloTree(*phyloTreeTest, phyloTreeName);
            maxLevel = phyloTree->getMaxLevel();
			
			m->mothurOut("DONE.\n"); 
			
			genusNodes = phyloTree->getGenusNodes(); 
			genusTotals = phyloTree->getGenusTotals();
			
            m->mothurOut("Reading template probabilities...     "); cout.flush();
            readProbFile(*probFileTest, *probFileTest2, probFileName, probFileName2);
			
        }else{
		
			//create search database and names vector
			generateDatabaseAndNames(txfile, tempFile, method, ksize, 0.0, 0.0, 0.0, 0.0, version);
			
			//prevents errors caused by creating shortcut files if you had an error in the sanity check.
			if (m->getControl_pressed()) {  util.mothurRemove(phyloTreeName);  util.mothurRemove(probFileName); util.mothurRemove(probFileName2); }
			else{ 
				genusNodes = phyloTree->getGenusNodes(); 
				genusTotals = phyloTree->getGenusTotals();
				
				m->mothurOut("Calculating template taxonomy tree...     "); cout.flush();
				
				phyloTree->printTreeNodes(phyloTreeName);
							
				m->mothurOut("DONE.\n"); 
				
				m->mothurOut("Calculating template probabilities...     "); cout.flush();
				
				numKmers = database->getMaxKmer() + 1;
			
				//initialze probabilities
				wordGenusProb.resize(numKmers);
                for (int j = 0; j < numKmers; j++) {  diffPair tempDiffPair; WordPairDiffArr.push_back(tempDiffPair); }
			
                for (int j = 0; j < wordGenusProb.size(); j++) {	wordGenusProb[j].resize(genusNodes.size(), 0.0);		}
                ofstream out; ofstream out2;

                if (shortcuts) { 
                    util.openOutputFile(probFileName, out); 
				
                    //output mothur version
                    out << "#" << version << endl;
				
                    out << numKmers << endl;
				
                    util.openOutputFile(probFileName2, out2);
				
                    //output mothur version
                    out2 << "#" << version << endl;
                }

				//for each word
				for (int i = 0; i < numKmers; i++) {
                    //m->mothurOut("[DEBUG]: kmer = " + toString(i) + "\n");
                    
					if (m->getControl_pressed()) {  break; }

                    if (shortcuts) {  out << i << '\t'; }
					
					vector<int> seqsWithWordi = database->getSequencesWithKmer(i);
					
					//for each sequence with that word
                    vector<int> count; count.resize(genusNodes.size(), 0);
					for (int j = 0; j < seqsWithWordi.size(); j++) {
						int temp = phyloTree->getGenusIndex(names[seqsWithWordi[j]]);
						count[temp]++;  //increment count of seq in this genus who have this word
					}
					
					//probabilityInTemplate = (# of seqs with that word in template + 0.50) / (total number of seqs in template + 1);
					float probabilityInTemplate = (seqsWithWordi.size() + 0.50) / (float) (names.size() + 1);
					diffPair tempProb(log(probabilityInTemplate), 0.0);
					WordPairDiffArr[i] = tempProb;
						
					int numNotZero = 0;
					for (int k = 0; k < genusNodes.size(); k++) {
						//probabilityInThisTaxonomy = (# of seqs with that word in this taxonomy + probabilityInTemplate) / (total number of seqs in this taxonomy + 1);
						
						
						wordGenusProb[i][k] = log((count[k] + probabilityInTemplate) / (float) (genusTotals[k] + 1));  
									
						if (count[k] != 0) {
                            if (shortcuts) { out << k << '\t' << wordGenusProb[i][k] << '\t' ; }
							numNotZero++;
						}
					}
					
                    
                    if (shortcuts) {
                        out << endl;
                        out2 << probabilityInTemplate << '\t' << numNotZero << '\t' << log(probabilityInTemplate) << endl;
                    }
                    
				}
                if (shortcuts) { out.close(); out2.close();  }
				
				//read in new phylotree with less info. - its faster
				ifstream phyloTreeTest(phyloTreeName.c_str());
				delete phyloTree;
				
				phyloTree = new PhyloTree(phyloTreeTest, phyloTreeName);
                maxLevel = phyloTree->getMaxLevel();
			}
		}
		
        if (m->getDebug()) { m->mothurOut("[DEBUG]: about to generateWordPairDiffArr\n"); }
		generateWordPairDiffArr();
        if (m->getDebug()) { m->mothurOut("[DEBUG]: done generateWordPairDiffArr\n"); }
        
        for (int i = 0; i < files.size(); i++) { delete files[i]; }
			
		m->mothurOut("DONE.\n");
		m->mothurOut("It took " + toString(time(nullptr) - start) + " seconds get probabilities.\n"); 
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "Bayesian");
		exit(1);
	}
}
/**************************************************************************************************/
Bayesian::~Bayesian() {
	try {
        if (phyloTree != nullptr) { delete phyloTree; }
        if (database != nullptr) {  delete database; }
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "~Bayesian");
		exit(1);
	}
}

/**************************************************************************************************/
string Bayesian::getTaxonomy(Sequence* seq, string& simpleTax, bool& flipped) {
	try {
		string tax = "";
        simpleTax = "";
		Kmer kmer(kmerSize);
		flipped = false;
		
		//get words contained in query
		//getKmerString returns a string where the index in the string is hte kmer number 
		//and the character at that index can be converted to be the number of times that kmer was seen
		string queryKmerString = kmer.getKmerString(seq->getUnaligned()); 
		
		vector<int> queryKmers;
		for (int i = 0; i < queryKmerString.length()-1; i++) {	// the -1 is to ignore any kmer with an N in it
			if (queryKmerString[i] != '!') { //this kmer is in the query
				queryKmers.push_back(i);
			}
		}
		
		//if user wants to test reverse compliment and its reversed use that instead
		if (flip) {	
			if (isReversed(queryKmers)) { 
				flipped = true;
				seq->reverseComplement(); 
				queryKmerString = kmer.getKmerString(seq->getUnaligned()); 
				queryKmers.clear();
				for (int i = 0; i < queryKmerString.length()-1; i++) {	// the -1 is to ignore any kmer with an N in it
					if (queryKmerString[i] != '!') { //this kmer is in the query
						queryKmers.push_back(i);
					}
				}
			}  
		}
		
		if (queryKmers.size() == 0) {  m->mothurOut(seq->getName() + " is bad. It has no kmers of length " + toString(kmerSize) + ".\n");  simpleTax = "unknown;";  return "unknown;"; }
		
		
		int index = getMostProbableTaxonomy(queryKmers);
		
		if (m->getControl_pressed()) { return tax; }
					
		//bootstrap - to set confidenceScore
        int numToSelect = queryKmers.size() / 8;
	
        if (m->getDebug()) {  m->mothurOut(seq->getName() + "\t"); }
        
		tax = bootstrapResults(queryKmers, index, numToSelect, simpleTax);
        
        if (m->getDebug()) {  m->mothurOut("\n"); }
		
		return tax;	
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "getTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/
string Bayesian::bootstrapResults(vector<int> kmers, int tax, int numToSelect, string& simpleTax) {
	try {
				
		map<int, int> confidenceScores; 
		
		//initialize confidences to 0 
		int seqIndex = tax;
		TaxNode seq = phyloTree->get(tax);
		confidenceScores[tax] = 0;
		
		while (seq.level != 0) { //while you are not at the root
			seqIndex = seq.parent;
			confidenceScores[seqIndex] = 0;
			seq = phyloTree->get(seq.parent);
		}
				
		map<int, int>::iterator itBoot;
		map<int, int>::iterator itBoot2;
		map<int, int>::iterator itConvert;
        
        int numKmers = kmers.size()-1;
        Utils util;
		for (int i = 0; i < iters; i++) {
			if (m->getControl_pressed()) { return "control"; }
			
            vector<int> temp;
			for (int j = 0; j < numToSelect; j++) {
				int index = util.getRandomIndex(numKmers);
                
				//add word to temp
				temp.push_back(kmers[index]);
			}
            
			//get taxonomy
			int newTax = getMostProbableTaxonomy(temp);
			//int newTax = 1;
			TaxNode taxonomyTemp = phyloTree->get(newTax);
			
			//add to confidence results
			while (taxonomyTemp.level != 0) { //while you are not at the root
				itBoot2 = confidenceScores.find(newTax); //is this a classification we already have a count on
				
				if (itBoot2 != confidenceScores.end()) { //this is a classification we need a confidence for
					(itBoot2->second)++;
				}
				
				newTax = taxonomyTemp.parent;
				taxonomyTemp = phyloTree->get(newTax);
			}
	
		}
		
		string confidenceTax = "";
		simpleTax = "";
		
		int seqTaxIndex = tax;
		TaxNode seqTax = phyloTree->get(tax);
		
        
		while (seqTax.level != 0) { //while you are not at the root
					
				itBoot2 = confidenceScores.find(seqTaxIndex); //is this a classification we already have a count on
				
				int confidence = 0;
				if (itBoot2 != confidenceScores.end()) { //already in confidence scores
					confidence = itBoot2->second;
				}
				
                if (m->getDebug()) { m->mothurOut(seqTax.name + "(" + toString(((confidence/(float)iters) * 100)) + ");"); }
            
				if (((confidence/(float)iters) * 100) >= confidenceThreshold) {
					confidenceTax = seqTax.name + "(" + toString(((confidence/(float)iters) * 100)) + ");" + confidenceTax;
					simpleTax = seqTax.name + ";" + simpleTax;
				}
            
				seqTaxIndex = seqTax.parent;
				seqTax = phyloTree->get(seqTax.parent);
		}
		
		if (confidenceTax == "") { confidenceTax = "unknown;"; simpleTax = "unknown;";  }
	
		return confidenceTax;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "bootstrapResults");
		exit(1);
	}
}
/**************************************************************************************************/
int Bayesian::getMostProbableTaxonomy(vector<int> queryKmer) {
	try {
		int indexofGenus = 0;
		
		double maxProbability = -1000000.0;
		//find taxonomy with highest probability that this sequence is from it

        for (int k = 0; k < genusNodes.size(); k++) {
			//for each taxonomy calc its probability
			
			double prob = 0.0000;
			for (int i = 0; i < queryKmer.size(); i++) { prob += wordGenusProb[queryKmer[i]][k]; }

			//is this the taxonomy with the greatest probability?
			if (prob > maxProbability) { 
				indexofGenus = genusNodes[k];
				maxProbability = prob;
			}
		}
			
		return indexofGenus;
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "getMostProbableTaxonomy");
		exit(1);
	}
}
//********************************************************************************************************************
//if it is more probable that the reverse compliment kmers are in the template, then we assume the sequence is reversed.
bool Bayesian::isReversed(vector<int>& queryKmers){
	try{
		bool reversed = false;
		float prob = 0;
		float reverseProb = 0;
		 
        for (int i = 0; i < queryKmers.size(); i++){
            int kmer = queryKmers[i];
            if (kmer >= 0){
                prob += WordPairDiffArr[kmer].prob;
				reverseProb += WordPairDiffArr[kmer].reverseProb;
            }
        }
		
        if (reverseProb > prob){ reversed = true; }
	
		return reversed;
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "isReversed");
		exit(1);
	}
}
//********************************************************************************************************************
int Bayesian::generateWordPairDiffArr(){
	try{
		Kmer kmer(kmerSize);
		for (int i = 0; i < WordPairDiffArr.size(); i++) {
			int reversedWord = kmer.getReverseKmerNumber(i);
			WordPairDiffArr[i].reverseProb = WordPairDiffArr[reversedWord].prob;
		}
		
		return 0;
	}catch(exception& e) {
		m->errorOut(e, "Bayesian", "generateWordPairDiffArr");
		exit(1);
	}
}
/**************************************************************************************************/
void Bayesian::readProbFile(ifstream& in, ifstream& inNum, string inName, string inNumName) {
	try{
		Utils util;
        //read version
        string line = util.getline(in); gobble(in);
        
        in >> numKmers; gobble(in);
        //initialze probabilities
        
        wordGenusProb.resize(numKmers);
        for (int j = 0; j < wordGenusProb.size(); j++) {	wordGenusProb[j].resize(genusNodes.size());		}
        
        int kmer, name, count;  count = 0;
        vector<int> num; num.resize(numKmers);
        float prob;
        vector<float> zeroCountProb; zeroCountProb.resize(numKmers);
        for (int j = 0; j < numKmers; j++) {  diffPair tempDiffPair; WordPairDiffArr.push_back(tempDiffPair); }
        
        //read version
        string line2 = util.getline(inNum); gobble(inNum);
        float probTemp;
        
        while (inNum) {
            inNum >> zeroCountProb[count] >> num[count] >> probTemp;
            WordPairDiffArr[count].prob = probTemp;
            count++;
            gobble(inNum);
            
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(zeroCountProb[count]) + '\t' + toString(num[count]) + '\t' + toString(numKmers) + "\n"); }

            
        }
        inNum.close();
 
        while(in) {
            in >> kmer;
            
            //set them all to zero value
            for (int i = 0; i < genusNodes.size(); i++) {
                wordGenusProb[kmer][i] = log(zeroCountProb[kmer] / (float) (genusTotals[i]+1));
            }
           
            //get probs for nonzero values
            for (int i = 0; i < num[kmer]; i++) {
                in >> name >> prob;
                wordGenusProb[kmer][name] = prob;
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + toString(name) + '\t' + toString(prob) + '\t' + toString(kmer) + "\n"); }
            }
            
            gobble(in);
        }
        in.close();
  		
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "readProbFile");
		exit(1);
	}
}
/**************************************************************************************************/






