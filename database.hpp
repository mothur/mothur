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

#include "mothur.h"

class Sequence;

class Database {
public:
	Database(string);
	virtual Sequence* findClosestSequence(Sequence*) = 0;
	virtual float getSearchScore();
protected:
	int numSeqs;
	float searchScore;
	vector<Sequence*> templateSequences;
};

#endif
