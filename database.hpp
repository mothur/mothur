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

/**************************************************************************************************/
struct seqMatch {  //used to select top n matches
		int seq;
		int match;
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
	Database();
	virtual ~Database();
	virtual void generateDB() = 0; 
	virtual void addSequence(Sequence) = 0;  //add sequence to search engine
	virtual vector<int> findClosestSequences(Sequence*, int) = 0;  // returns indexes of n closest sequences to query
	virtual vector<int> findClosestMegaBlast(Sequence*, int){return results;}
	virtual float getSearchScore();
	virtual int getLongestBase(); 
	virtual void readKmerDB(ifstream&){};
	virtual void setNumSeqs(int i) {	numSeqs = i; 	}
	virtual vector<int> getSequencesWithKmer(int){ vector<int> filler; return filler; };  
	virtual int getMaxKmer(){	return 1;	};

	
protected:
	MothurOut* m;
	int numSeqs, longest;
	float searchScore;
	vector<int> results;
};
/**************************************************************************************************/
#endif
