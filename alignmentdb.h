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

	AlignmentDB(string, string, int, int, int, int, int);  //reads fastafile passed in and stores sequences
	~AlignmentDB();
	
	Sequence findClosestSequence(Sequence*);
	float getSearchScore()  {  return search->getSearchScore();  }
	int getLongestBase()	{  return longest;  }
	
private:
	int numSeqs, longest;
	float searchScore;
	
	Database* search;
	vector<Sequence> templateSequences;
	Sequence emptySequence;
};

/**************************************************************************************************/

#endif

