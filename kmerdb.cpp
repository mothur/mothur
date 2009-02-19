/*
 *  kmerdb.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;

#include "mothur.h"
#include "sequence.hpp"
#include "kmer.hpp"
#include "database.hpp"
#include "kmerdb.hpp"

/**************************************************************************************************/

KmerDB::KmerDB(string fastaFileName, int kSize) : Database(fastaFileName), kmerSize(kSize) {

	string kmerDBName = fastaFileName.substr(0,fastaFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
	ifstream kmerFileTest(kmerDBName.c_str());
	
	int power4s[9] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536 };
	
	maxKmer = power4s[kmerSize];
	kmerLocations.resize(maxKmer+1);
	
	if(!kmerFileTest){
		cout << "Generating the " << kmerDBName << " database...\t";	cout.flush();
		generateKmerDB(kmerDBName);	
	}
	else{
		cout << "Reading in the " << kmerDBName << " database...\t";	cout.flush();
		readKmerDB(kmerDBName, kmerFileTest);
	}
	cout << "DONE." << endl << endl;	cout.flush();

}

/**************************************************************************************************/

Sequence* KmerDB::findClosestSequence(Sequence* candidateSeq){
	
	vector<int> matches(numSeqs, 0);
	vector<int> timesKmerFound(kmerLocations.size()+1, 0);
	
	int maxMatches = 0;
	int maxSequence = 0;
	
	string query = candidateSeq->getUnaligned();
	
	int numKmers = query.length() - kmerSize + 1;
	Kmer kmer(kmerSize);
	
	for(int i=0;i<numKmers;i++){
		
		int kmerNumber = kmer.getKmerNumber(query, i);
		
		if(timesKmerFound[kmerNumber] == 0){
			for(int j=0;j<kmerLocations[kmerNumber].size();j++){
				matches[kmerLocations[kmerNumber][j]]++;
			}
		}
		timesKmerFound[kmerNumber] = 1;
		
	}
	for(int i=0;i<numSeqs;i++){
		if(matches[i] > maxMatches){
			maxMatches = matches[i];
			maxSequence = i;
		}
	}
	return templateSequences[maxSequence];
	
}

/**************************************************************************************************/

void KmerDB::generateKmerDB(string kmerDBName){
	
	
	Kmer kmer(kmerSize);
	
	for(int i=0;i<numSeqs;i++){

		string seq = templateSequences[i]->getUnaligned();
		int numKmers = seq.length() - kmerSize + 1;
		
		for(int j=0;j<numKmers;j++){
			int kmerNumber = kmer.getKmerNumber(seq, j);
			kmerLocations[kmerNumber].push_back(i);
		}
	}
	
	ofstream kmerFile(kmerDBName.c_str(), ios::trunc);
	if(!kmerFile) {
		cerr << "Error: Could not open " << kmerDBName << endl;
		exit(1);
	}
	
	for(int i=0;i<maxKmer;i++){
		kmerFile << i << ' ' << kmerLocations[i].size();
		for(int j=0;j<kmerLocations[i].size();j++){
			kmerFile << ' ' << kmerLocations[i][j];
		}
		kmerFile << endl;
	}
	kmerFile.close();
	
}

/**************************************************************************************************/

void KmerDB::readKmerDB(string kmerDBName, ifstream& kmerDBFile){

	kmerDBFile.seekg(0);
	
	string seqName;
	int seqNumber;
	
	for(int i=0;i<numSeqs;i++){
		int numValues;
		kmerDBFile >> seqName >> numValues;
		
		for(int j=0;j<numValues;j++){
			kmerDBFile >> seqNumber;
			kmerLocations[i].push_back(seqNumber);
		}
	}
	kmerDBFile.close();
}

/**************************************************************************************************/
