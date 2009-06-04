#ifndef DISTANCEDB_HPP
#define DISTANCEDB_HPP

/*
 *  distancedb.hpp
 *  
 *
 *  Created by Pat Schloss on 12/29/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


#include "mothur.h"

class DistanceDB : public Database {
	
public:
	DistanceDB(string, string);
	Sequence* findClosestSequence(Sequence*);
	
private:

	struct hit{
		string seqName;
		int indexNumber;
		float simScore;
	};
	vector<hit> mostSimSequenceVector;
//	ifstream inputData;
	int searchIndex;
};

#endif
