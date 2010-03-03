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

/**************************************************************************************************/
Bayesian::Bayesian(string tfile, string tempFile, string method, int ksize, int cutoff, int i) : 
Classify(tfile, tempFile, method, ksize, 0.0, 0.0, 0.0, 0.0), kmerSize(ksize), confidenceThreshold(cutoff), iters(i)  {
	try {
					
		numKmers = database->getMaxKmer() + 1;
		
		//initialze probabilities
		wordGenusProb.resize(numKmers);
		
		genusNodes = phyloTree->getGenusNodes(); 
		
		for (int j = 0; j < wordGenusProb.size(); j++) {	wordGenusProb[j].resize(genusNodes.size());		}
			
		//reset counts because we are on a new word
		for (int j = 0; j < genusNodes.size(); j++) {
			TaxNode temp = phyloTree->get(genusNodes[j]);
			genusTotals.push_back(temp.accessions.size());
		}

		
		/************calculate the probablity that each word will be in a specific taxonomy*************/
		ofstream out;
		string probFileName = tempFile.substr(0,tempFile.find_last_of(".")+1) + char('0'+ kmerSize) + "mer.prob";
		ifstream probFileTest(probFileName.c_str());
		
		ofstream out2;
		string probFileName2 = tempFile.substr(0,tempFile.find_last_of(".")+1) + char('0'+ kmerSize) + "mer.numNonZero";
		ifstream probFileTest2(probFileName2.c_str());
		
		int start = time(NULL);
		
		if(probFileTest && probFileTest2){	
			m->mothurOut("Reading template probabilities...     "); cout.flush();
			readProbFile(probFileTest, probFileTest2);	
		}else{
			m->mothurOut("Calculating template probabilities...     "); cout.flush();

			ofstream out;
			openOutputFile(probFileName, out);
			
			ofstream out2;
			openOutputFile(probFileName2, out2);
			
			//for each word
			for (int i = 0; i < numKmers; i++) {
				if (m->control_pressed) { break; }
				
				out << i << '\t';
				
				vector<int> seqsWithWordi = database->getSequencesWithKmer(i);
				
				map<int, int> count;
				for (int k = 0; k < genusNodes.size(); k++) {  count[genusNodes[k]] = 0;  }			
						
				//for each sequence with that word
				for (int j = 0; j < seqsWithWordi.size(); j++) {
					int temp = phyloTree->getIndex(names[seqsWithWordi[j]]);
					count[temp]++;  //increment count of seq in this genus who have this word
				}
				
				//probabilityInTemplate = (# of seqs with that word in template + 0.05) / (total number of seqs in template + 1);
				float probabilityInTemplate = (seqsWithWordi.size() + 0.50) / (float) (names.size() + 1);
				
				int numNotZero = 0;
				for (int k = 0; k < genusNodes.size(); k++) {
					//probabilityInThisTaxonomy = (# of seqs with that word in this taxonomy + probabilityInTemplate) / (total number of seqs in this taxonomy + 1);
					wordGenusProb[i][k] = log((count[genusNodes[k]] + probabilityInTemplate) / (float) (genusTotals[k] + 1));  
					if (count[genusNodes[k]] != 0) {  out << k << '\t' << wordGenusProb[i][k] << '\t';  numNotZero++;  }
				}
				out << endl;
				out2 << probabilityInTemplate << '\t' << numNotZero << endl;
			}
			
			out.close();
			out2.close();
		}
		
		
		m->mothurOut("DONE."); m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " seconds get probabilities. "); m->mothurOutEndLine();
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "Bayesian");
		exit(1);
	}
}
/**************************************************************************************************/
string Bayesian::getTaxonomy(Sequence* seq) {
	try {
		string tax = "";
		Kmer kmer(kmerSize);
		
		//get words contained in query
		//getKmerString returns a string where the index in the string is hte kmer number 
		//and the character at that index can be converted to be the number of times that kmer was seen
		string queryKmerString = kmer.getKmerString(seq->getUnaligned()); 
		vector<int> queryKmers;
		for (int i = 0; i < queryKmerString.length(); i++) {
			if (queryKmerString[i] != '!') { //this kmer is in the query
				queryKmers.push_back(i);
			}
		}
	
		int index = getMostProbableTaxonomy(queryKmers);
		
		if (m->control_pressed) { return tax; }
					
		//bootstrap - to set confidenceScore
		int numToSelect = queryKmers.size() / 8;
		tax = bootstrapResults(queryKmers, index, numToSelect);
						
		return tax;	
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "getTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/
string Bayesian::bootstrapResults(vector<int> kmers, int tax, int numToSelect) {
	try {
		
		//taxConfidenceScore.clear(); //clear out previous seqs scores
				
		vector< map<string, int> > confidenceScores; //you need the added vector level of confusion to account for the level that name is seen since they can be the same
										//map of classification to confidence for all areas seen
									   //ie. Bacteria;Alphaproteobacteria;Rhizobiales;Azorhizobium_et_rel.;Methylobacterium_et_rel.;Bosea;
									   //ie. Bacteria -> 100, Alphaproteobacteria -> 100, Rhizobiales -> 87, Azorhizobium_et_rel. -> 78, Methylobacterium_et_rel. -> 70, Bosea -> 50
		confidenceScores.resize(100);  //if you have more than 100 levels of classification...
		
		map<string, int>::iterator itBoot;
		map<string, int>::iterator itBoot2;
		map<int, int>::iterator itConvert;
		
		for (int i = 0; i < iters; i++) {
			if (m->control_pressed) { return "control"; }
			
			vector<int> temp;
						
			for (int j = 0; j < numToSelect; j++) {
				int index = int(rand() % kmers.size());
				
				//add word to temp
				temp.push_back(kmers[index]);
			}
			
			//get taxonomy
			int newTax = getMostProbableTaxonomy(temp);
			TaxNode taxonomy = phyloTree->get(newTax);
			
			//add to confidence results
			while (taxonomy.level != 0) { //while you are not at the root
				
				itBoot2 = confidenceScores[taxonomy.level].find(taxonomy.name); //is this a classification we already have a count on
				
				if (itBoot2 == confidenceScores[taxonomy.level].end()) { //not already in confidence scores
					confidenceScores[taxonomy.level][taxonomy.name] = 1;
				}else{
					confidenceScores[taxonomy.level][taxonomy.name]++;
				}
			
				taxonomy = phyloTree->get(taxonomy.parent);
			}
		}
		
		string confidenceTax = "";
		simpleTax = "";
		TaxNode seqTax = phyloTree->get(tax);
		
		while (seqTax.level != 0) { //while you are not at the root
				
				itBoot2 = confidenceScores[seqTax.level].find(seqTax.name); //is this a classification we already have a count on
				
				int confidence = 0;
				if (itBoot2 != confidenceScores[seqTax.level].end()) { //not already in confidence scores
					confidence = confidenceScores[seqTax.level][seqTax.name];
				}
				
				if (confidence >= confidenceThreshold) {
					confidenceTax = seqTax.name + "(" + toString(((confidence/(float)iters) * 100)) + ");" + confidenceTax;
					simpleTax = seqTax.name + ";" + simpleTax;
				}
				
				seqTax = phyloTree->get(seqTax.parent);
		}
		
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
		int indexofGenus;
		
		double maxProbability = -1000000.0;
		//find taxonomy with highest probability that this sequence is from it
		for (int k = 0; k < genusNodes.size(); k++) {
		
			//for each taxonomy calc its probability
			double prob = 1.0;
			for (int i = 0; i < queryKmer.size(); i++) {
				prob += wordGenusProb[queryKmer[i]][k];
			}
		
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
/*************************************************************************************************
map<string, int> Bayesian::parseTaxMap(string newTax) {
	try{
	
		map<string, int> parsed;
		
		newTax = newTax.substr(0, newTax.length()-1);  //get rid of last ';'
	
		//parse taxonomy
		string individual;
		while (newTax.find_first_of(';') != -1) {
			individual = newTax.substr(0,newTax.find_first_of(';'));
			newTax = newTax.substr(newTax.find_first_of(';')+1, newTax.length());
			parsed[individual] = 1;
		}
		
		//get last one
		parsed[newTax] = 1;

		return parsed;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Bayesian", "parseTax");
		exit(1);
	}
}
/**************************************************************************************************/
void Bayesian::readProbFile(ifstream& in, ifstream& inNum) {
	try{
		
		int kmer, name, count;  count = 0;
		vector<int> num; num.resize(numKmers);
		float prob;
		vector<float> zeroCountProb; zeroCountProb.resize(numKmers);		
		
		while (inNum) {
			inNum >> zeroCountProb[count] >> num[count];  
			count++;
			gobble(inNum);
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






