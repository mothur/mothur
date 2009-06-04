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

KmerDB::KmerDB(string fastaFileName, int kSize) : Database(fastaFileName), kmerSize(kSize) {

	string kmerDBName = fastaFileName.substr(0,fastaFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
	ifstream kmerFileTest(kmerDBName.c_str());
	
	int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
	
	maxKmer = power4s[kmerSize];
	kmerLocations.resize(maxKmer+1);
	
	if(!kmerFileTest){		//	if we can open the kmer db file, then read it in...
		cout << "Generating the " << kmerDBName << " database...\t";	cout.flush();
		generateKmerDB(kmerDBName);	
	}
	else{					//	...otherwise generate it.
		cout << "Reading in the " << kmerDBName << " database...\t";	cout.flush();
		readKmerDB(kmerDBName, kmerFileTest);
	}
	cout << "DONE." << endl << endl;	cout.flush();

}

/**************************************************************************************************/

Sequence* KmerDB::findClosestSequence(Sequence* candidateSeq){
	
	Kmer kmer(kmerSize);
	
	searchScore = 0;
	int maxSequence = 0;

	vector<int> matches(numSeqs, 0);						//	a record of the sequences with shared kmers
	vector<int> timesKmerFound(kmerLocations.size()+1, 0);	//	a record of the kmers that we have already found

	int numKmers = candidateSeq->getNumBases() - kmerSize + 1;	

	for(int i=0;i<numKmers;i++){
		int kmerNumber = kmer.getKmerNumber(candidateSeq->getUnaligned(), i);		//	go through the query sequence and get a kmer number
		if(timesKmerFound[kmerNumber] == 0){				//	if we haven't seen it before...
			for(int j=0;j<kmerLocations[kmerNumber].size();j++){//increase the count for each sequence that also has
				matches[kmerLocations[kmerNumber][j]]++;	//	that kmer
			}
		}
		timesKmerFound[kmerNumber] = 1;						//	ok, we've seen the kmer now
	}

	for(int i=0;i<numSeqs;i++){								//	find the db sequence that shares the most kmers with
		if(matches[i] > searchScore){					//	the query sequence
			searchScore = matches[i];
			maxSequence = i;
		}
	}

	searchScore = 100 * searchScore / (float) numKmers;		//	return the Sequence object corresponding to the db
	return templateSequences[maxSequence];					//	sequence with the most shared kmers.
}

/**************************************************************************************************/

void KmerDB::generateKmerDB(string kmerDBName){
	
	Kmer kmer(kmerSize);
	
	for(int i=0;i<numSeqs;i++){								//	for all of the template sequences...

		string seq = templateSequences[i]->getUnaligned();	//	...take the unaligned sequence...
		int numKmers = seq.length() - kmerSize + 1;
		
		vector<int> seenBefore(maxKmer+1,0);
		for(int j=0;j<numKmers;j++){						//	...step though the sequence and get each kmer...
			int kmerNumber = kmer.getKmerNumber(seq, j);
			if(seenBefore[kmerNumber] == 0){
				kmerLocations[kmerNumber].push_back(i);		//	...insert the sequence index into kmerLocations for
			}												//	the appropriate kmer number
			seenBefore[kmerNumber] = 1;
		}													
	}
	
	ofstream kmerFile;										//	once we have the kmerLocations folder print it out
	openOutputFile(kmerDBName, kmerFile);					//	to a file
	
	for(int i=0;i<maxKmer;i++){								//	step through all of the possible kmer numbers
		kmerFile << i << ' ' << kmerLocations[i].size();	//	print the kmer number and the number of sequences with
		for(int j=0;j<kmerLocations[i].size();j++){			//	that kmer.  then print out the indices of the sequences
			kmerFile << ' ' << kmerLocations[i][j];			//	with that kmer.
		}
		kmerFile << endl;
	}
	kmerFile.close();
	
}

/**************************************************************************************************/

void KmerDB::readKmerDB(string kmerDBName, ifstream& kmerDBFile){

	kmerDBFile.seekg(0);									//	start at the beginning of the file
	
	string seqName;
	int seqNumber;
	
	for(int i=0;i<maxKmer;i++){
		int numValues;	
		kmerDBFile >> seqName >> numValues;
		
		for(int j=0;j<numValues;j++){						//	for each kmer number get the...
			kmerDBFile >> seqNumber;						//		1. number of sequences with the kmer number
			kmerLocations[i].push_back(seqNumber);			//		2. sequence indices
		}
	}
	kmerDBFile.close();
}

/**************************************************************************************************/
