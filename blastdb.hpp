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
	BlastDB(float, float, float, float);
	BlastDB();
	~BlastDB();
	
	void generateDB();
	void addSequence(Sequence);
	vector<int> findClosestSequences(Sequence*, int);
	vector<int> findClosestMegaBlast(Sequence*, int, int);
	
private:
	string dbFileName;
	string queryFileName;
	string blastFileName;
	string path;
	
	int count;
	float gapOpen;
	float gapExtend;
	float match;
	float misMatch;
	
};

#endif
