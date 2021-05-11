/*
 *  sequencedb.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sequencedb.h"
#include "sequence.hpp"
#include "mothur.h"
#include "calculator.h"
#include "kmer.hpp"

/***********************************************************************/

SequenceDB::SequenceDB() {  m = MothurOut::getInstance();  length = 0; samelength = true; }
/***********************************************************************/
//the clear function free's the memory
SequenceDB::~SequenceDB() { data.clear(); }

/***********************************************************************/

SequenceDB::SequenceDB(int newSize) {
	data.resize(newSize, Sequence());
	length = 0; samelength = true;
}

/***********************************************************************/
//kmerDB[0] = vector<int> maxKmers long, contains kmer counts
SequenceDB::SequenceDB(ifstream& filehandle, int kmerSize, vector< vector< int > >& kmerDB, vector< int >& lengths) {
    try{
        Utils util; length = 0; samelength = true; lengths.clear(); kmerDB.clear();
        
        int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };
        
        int maxKmer = power4s[kmerSize];
        
        Kmer kmer(kmerSize);
        
        while (!filehandle.eof()) {
            //input sequence info into sequencedb
            Sequence newSequence(filehandle); util.gobble(filehandle);
            
            if (newSequence.getName() != "") {
                if (length == 0) { length = newSequence.getAligned().length(); }
                if (length != newSequence.getAligned().length()) { samelength = false;  }
                
                data.push_back(newSequence);
                
                vector<int> kmerLocations; kmerLocations.resize(maxKmer+1, 0);
                
                int numKmers = newSequence.getNumBases() - kmerSize + 1;
            
                for(int i=0;i<numKmers;i++){
                    int kmerNumber = kmer.getKmerNumber(newSequence.getUnaligned(), i);
                    
                    kmerLocations[kmerNumber]++; //ok, we've seen the kmer now
                }
                
                kmerDB.push_back(kmerLocations);
                lengths.push_back(newSequence.getNumBases());
            }
        }

        filehandle.close();
        
    }
    catch(exception& e) {
        m->errorOut(e, "SequenceDB", "SequenceDB");
        exit(1);
    }
}
/***********************************************************************/
SequenceDB::SequenceDB(ifstream& filehandle) {
	try{
		length = 0; samelength = true;
        Utils util;
		//read through file
		while (!filehandle.eof()) {
			//input sequence info into sequencedb
			Sequence newSequence(filehandle);
			
			if (newSequence.getName() != "") {   
				if (length == 0) { length = newSequence.getAligned().length(); }
                if (length != newSequence.getAligned().length()) { samelength = false;  }
				data.push_back(newSequence);  
			}
			
			//takes care of white space
			util.gobble(filehandle);
		}

		filehandle.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "SequenceDB");
		exit(1);
	}
}
/***********************************************************************/

int SequenceDB::getNumSeqs() {
	return data.size();
}

/***********************************************************************/
Sequence SequenceDB::get(int index) {
	return data[index];
}
/***********************************************************************/

void SequenceDB::push_back(Sequence newSequence) {
	try {
		if (length == 0) { length = newSequence.getAligned().length(); }
		if (length != newSequence.getAligned().length()) { samelength = false; }

		data.push_back(newSequence);
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "push_back");
		exit(1);
	}
}

/***********************************************************************/

