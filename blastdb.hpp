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
	BlastDB(string, float, float, float, float, string);
	BlastDB(string);
	BlastDB(const BlastDB& bdb) : dbFileName(bdb.dbFileName), queryFileName(bdb.queryFileName), blastFileName(bdb.blastFileName), path(bdb.path),
									count(bdb.count), gapOpen(bdb.gapOpen), gapExtend(bdb.gapExtend), match(bdb.match), misMatch(bdb.misMatch), Database(bdb) {}
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
