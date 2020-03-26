#ifndef DATABASE_HPP
#define DATABASE_HPP

/*
 *  database.hpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


/* This class is a parent to blastdb, kmerdb, suffixdb.  */

#include "mothur.h"
#include "sequence.hpp"
#include "currentfile.h"
#include "utils.hpp"

/**************************************************************************************************/
struct classifierOTU {
    vector<vector<char> > otuData; //otuData[0] -> vector of first characters from each sequence in the OTU, otuData[1] -> vector of second characters from each sequence in the OTU
    
    /*
     otuData.size = num columns in seq's alignment
     otuData[i].size() = numSeqs in otu
     
     seq1 > atgcaag
     seq2 > gacctga
     seq3 > cctgacg
     
     otuData[0] > {a,g,c}
     otuData[1] > {t,a,c}
     otuData[2] > {g,c,t}
     
     otuData[i] > {charInAllCols} if all chars in otuData[i] are identical. ie, ignore column
     otuData[i] > {a} all seqs contain 'a' in column i of alignment
     */
    
    vector<int> bootstrapColumns; // initialize to all columns
    
    classifierOTU() {}
    classifierOTU(vector<vector<char> > otu) : otuData(otu) {
        for (int i = 0; i < otu.size(); i++) { bootstrapColumns.push_back(i); }
    }
};
/**************************************************************************************************/
struct seqMatch {  //used to select top n matches
		int seq;
		int match;
		seqMatch() {}
		seqMatch(int s, int m) : seq(s), match(m) {}
};
/**************************************************************************************************/
inline bool compareSeqMatches (seqMatch member, seqMatch member2){ //sorts largest to smallest
	if(member.match > member2.match){
		return true;   }   
	else{
		return false; 
	}
}
/**************************************************************************************************/
inline bool compareSeqMatchesReverse (seqMatch member, seqMatch member2){ //sorts largest to smallest
	if(member.match < member2.match){
		return true;   }   
	else{
		return false; 
	}
}

/**************************************************************************************************/
class Database {

public:
    Database(){ longest = 0; numSeqs = 0; m = MothurOut::getInstance();  }
    
    virtual ~Database(){};
    virtual void generateDB() = 0;
    virtual void readKmerDB(ifstream&){};
    virtual void addSequence(Sequence) = 0;  //add sequence to search engine
    virtual void addSequences(vector<Sequence> seqs) { for (int i = 0; i < seqs.size(); i++) { addSequence(seqs[i]); } }

	virtual void setNumSeqs(int i) {	numSeqs = i; 	}
    
    
	virtual vector<int> findClosestSequences(Sequence*, int, vector<float>&) const = 0;  // returns indexes of n closest sequences to query
	
    virtual int getLongestBase() {	return longest+1;		}
    virtual vector<int> getSequencesWithKmer(int){ vector<int> filler; return filler; };
	virtual int getReversed(int) { return 0; } 
	virtual int getMaxKmer(){	return 1;	}
    virtual string getName(int) { return ""; }
    virtual vector<int> findClosestMegaBlast(Sequence*, int, int){ vector<int> results; return results;}
	
protected:
    
    MothurOut* m;
	int numSeqs, longest;
    Utils util;
	
	
};
/**************************************************************************************************/
#endif
