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

class Database {
public:
	Database(string);
	virtual Sequence* findClosestSequence(Sequence*) = 0;
	
protected:
	int numSeqs;
	vector<Sequence*> templateSequences;	
};

#endif
