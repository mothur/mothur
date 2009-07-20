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


/* This class is a parent to blastdb, distancedb, kmerdb, suffixdb.  Which are used to convert a squencedb object into that form. */

#include "mothur.h"

class Sequence;

class Database {
public:
	Database(string);
	virtual ~Database();
	virtual Sequence findClosestSequence(Sequence*) = 0;
	virtual float getSearchScore();
	virtual int getLongestBase(); 
	
protected:
	int numSeqs, longest;
	float searchScore;
	vector<Sequence> templateSequences;
};

#endif
