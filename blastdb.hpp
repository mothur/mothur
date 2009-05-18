#ifndef BLASTDB_HPP
#define BLASTDB_HPP


/*
 *  blastdb.hpp
 *  
 *
 *  Created by Pat Schloss on 12/22/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"

class BlastDB : public Database {
public:
	BlastDB(string, float, float, float, float);
	~BlastDB();
	Sequence* findClosestSequence(Sequence*);

private:
	string dbFileName;
	string queryFileName;
	string blastFileName;
	
	float gapOpen;
	float gapExtend;
	float match;
	float misMatch;
	Sequence* emptySequence;
};

#endif
