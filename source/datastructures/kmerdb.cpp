/*
 *  kmerdb.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This class is a child class of the Database class, which stores the template sequences as a kmer table and provides
 *	a method of searching the kmer table for the sequence with the most kmers in common with a query sequence.
 *	kmerLocations is the primary storage variable that is a two-dimensional vector where each row represents the
 *	different number of kmers and each column contains the index to sequences that use that kmer.
 *
 *	Construction of an object of this type will first look for an appropriately named database file and if it is found
 *	then will read in the database file (readKmerDB), otherwise it will generate one and store the data in memory
 *	(generateKmerDB)
 *
 *	The search method used here is roughly the same as that used in the SimRank program that is found at the
 *	greengenes website.  The default kmer size is 7.  The speed complexity is between O(L) and O(LN).  When I use 7mers
 *	on average a kmer is found in ~100 other sequences with a database of ~5000 sequences.  If this is the case then the
 *	time would be on the order of O(0.1LN) -> fast.
 *

 */

#include "sequence.hpp"
#include "kmer.hpp"
#include "database.hpp"
#include "kmerdb.hpp"

/**************************************************************************************************/

KmerDB::KmerDB(string fastaFileName, int kSize) : Database(), kmerSize(kSize) {
	try { 
	
		kmerDBName = fastaFileName.substr(0,fastaFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
		
		int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
		count = 0;
		
		maxKmer = power4s[kmerSize];
		kmerLocations.resize(maxKmer+1);
        
        CurrentFile* current; current = CurrentFile::getInstance();
        version = current->getVersion();
		
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "KmerDB");
		exit(1);
	}	

}
/**************************************************************************************************/
KmerDB::KmerDB() : Database() {
    CurrentFile* current; current = CurrentFile::getInstance();
    version = current->getVersion();
}
/**************************************************************************************************/

KmerDB::~KmerDB(){}

/**************************************************************************************************/

vector<int> KmerDB::findClosestSequences(Sequence* candidateSeq, int num, vector<float>& Scores) const{
	try {
		if (num > numSeqs) { m->mothurOut("[WARNING]: you requested " + toString(num) + " closest sequences, but the template only contains " + toString(numSeqs) + ", adjusting.\n");  num = numSeqs; }
		
		vector<int> topMatches;
		Kmer kmer(kmerSize);
		float searchScore = 0;
		Scores.clear();
		
		vector<int> matches(numSeqs, 0);						//	a record of the sequences with shared kmers
		vector<bool> timesKmerFound(kmerLocations.size()+1, false);	//	a record of the kmers that we have already found
		
		int numKmers = candidateSeq->getNumBases() - kmerSize + 1;	
	
		for(int i=0;i<numKmers;i++){
			int kmerNumber = kmer.getKmerNumber(candidateSeq->getUnaligned(), i);		//	go through the query sequence and get a kmer number
			if(!timesKmerFound[kmerNumber]){				//	if we haven't seen it before...
				for(int j=0;j<kmerLocations[kmerNumber].size();j++){//increase the count for each sequence that also has
					matches[kmerLocations[kmerNumber][j]]++;	//	that kmer
				}
			}
			timesKmerFound[kmerNumber] = true;						//	ok, we've seen the kmer now
		}
		
		if (num != 1) {
			vector<seqMatch> seqMatches; seqMatches.resize(numSeqs);
			for(int i=0;i<numSeqs;i++){		
				seqMatches[i].seq = i;
				seqMatches[i].match = matches[i];
			}
			
			//sorts putting largest matches first
			sort(seqMatches.begin(), seqMatches.end(), compareSeqMatches);
			
			searchScore = seqMatches[0].match;
			searchScore = 100 * searchScore / (float) numKmers;		//	return the Sequence object corresponding to the db
            Scores.push_back(searchScore);
            
			//save top matches
			for (int i = 0; i < num; i++) {
				topMatches.push_back(seqMatches[i].seq);
				float thisScore = 100 * seqMatches[i].match / (float) numKmers;
				Scores.push_back(thisScore);
			}
		}else{
			int bestIndex = 0;
			int bestMatch = -1;
            
			for(int i=0;i<numSeqs;i++){	
				
				if (matches[i] > bestMatch) {
					bestIndex = i;
					bestMatch = matches[i];
				}
			}
            
			searchScore = bestMatch;
			searchScore = 100 * searchScore / (float) numKmers;		//	return the Sequence object corresponding to the db
			topMatches.push_back(bestIndex);
			Scores.push_back(searchScore);
		}
		return topMatches;		
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "findClosestSequences");
		exit(1);
	}	
}

/**************************************************************************************************/
//print shortcut file
void KmerDB::generateDB(){
	try {
		
		ofstream kmerFile;										//	once we have the kmerLocations folder print it out
		util.openOutputFile(kmerDBName, kmerFile);					//	to a file
		
		//output version
		kmerFile << "#" << version << endl;
		
		for(int i=0;i<maxKmer;i++){								//	step through all of the possible kmer numbers
			kmerFile << i << ' ' << kmerLocations[i].size();	//	print the kmer number and the number of sequences with
			for(int j=0;j<kmerLocations[i].size();j++){			//	that kmer.  then print out the indices of the sequences
				kmerFile << ' ' << kmerLocations[i][j];			//	with that kmer.
			}
			kmerFile << endl;
		}
		kmerFile.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "generateDB");
		exit(1);
	}	
	
}
/**************************************************************************************************/
void KmerDB::addSequence(Sequence seq) {
	try {
		Kmer kmer(kmerSize);
		
		string unaligned = seq.getUnaligned();	//	...take the unaligned sequence...
		int numKmers = unaligned.length() - kmerSize + 1;
			
		vector<int> seenBefore(maxKmer+1,0);
		for(int j=0;j<numKmers;j++){						//	...step though the sequence and get each kmer...
			int kmerNumber = kmer.getKmerNumber(unaligned, j);
			if(seenBefore[kmerNumber] == 0){
				kmerLocations[kmerNumber].push_back(count);		//	...insert the sequence index into kmerLocations for
			}												//	the appropriate kmer number
			seenBefore[kmerNumber] = 1;
		}													
	
		count++;
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "addSequence");
		exit(1);
	}	
}
/**************************************************************************************************/

void KmerDB::readKmerDB(ifstream& kmerDBFile){
	try {
					
		kmerDBFile.seekg(0);									//	start at the beginning of the file
		
		//read version
		string line = util.getline(kmerDBFile); util.gobble(kmerDBFile);
		
		string seqName;
		int seqNumber;

		for(int i=0;i<maxKmer;i++){
			int numValues = 0;	
			kmerDBFile >> seqName >> numValues;
			
			for(int j=0;j<numValues;j++){						//	for each kmer number get the...
				kmerDBFile >> seqNumber;						//		1. number of sequences with the kmer number
				kmerLocations[i].push_back(seqNumber);			//		2. sequence indices
			}
		}
		kmerDBFile.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "readKmerDB");
		exit(1);
	}	
}

/**************************************************************************************************/
int KmerDB::getCount(int kmer) {
	try {
		if (kmer < 0) { return 0; }  //if user gives negative number
		else if (kmer > maxKmer) {	return 0;	}  //or a kmer that is bigger than maxkmer
		else {	return kmerLocations[kmer].size();	}  // kmer is in vector range
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "getCount");
		exit(1);
	}	
}
/**************************************************************************************************/
int KmerDB::getReversed(int kmerNumber) {
	try {
		Kmer kmer(kmerSize);
		
		if (kmerNumber < 0) { return 0; }  //if user gives negative number
		else if (kmerNumber > maxKmer) {	return 0;	}  //or a kmer that is bigger than maxkmer
		else {	return kmer.getReverseKmerNumber(kmerNumber);	}  // kmer is in vector range
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "getReversed");
		exit(1);
	}	
}
/**************************************************************************************************/
vector<int> KmerDB::getSequencesWithKmer(int kmer) {
	try {
		
		vector<int> seqs;
	
		if (kmer < 0) { }  //if user gives negative number
		else if (kmer > maxKmer) {	}  //or a kmer that is bigger than maxkmer
		else {	seqs = kmerLocations[kmer];	}
		
		return seqs;
	}
	catch(exception& e) {
		m->errorOut(e, "KmerDB", "getSequencesWithKmer");
		exit(1);
	}	
}
/**************************************************************************************************/


/**************************************************************************************************/


