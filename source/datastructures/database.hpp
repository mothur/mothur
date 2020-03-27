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
