#ifndef ALIGNMENTDB_H
#define ALIGNMENTDB_H

/*
 *  alignmentdb.h
 *  Mothur
 *
 *  Created by westcott on 11/4/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "sequence.hpp"
#include "database.hpp"

/**************************************************************************************************/

class AlignmentDB {

public:

	AlignmentDB(string, string, int, float, float, float, float, int, bool);  //reads fastafile passed in and stores sequences
	AlignmentDB(string);
	~AlignmentDB();
	
	Sequence findClosestSequence(Sequence*, float&); //sequence to align, searchScore
	int getLongestBase()	{  return longest;  }
	
private:
    std::mutex mutex;
	int numSeqs, longest, threadID;
	string method;
	
	Database* search;
	vector<Sequence> templateSequences;
	Sequence emptySequence;
	MothurOut* m;
};

/**************************************************************************************************/

#endif

