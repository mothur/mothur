#ifndef DISTANCEDB_HPP
#define DISTANCEDB_HPP

/*
 *  distancedb.hpp
 *  
 *
 *  Created by westcott on 1/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "calculator.h"

class DistanceDB : public SearchDatabase {
	
public:
	
	DistanceDB();
	~DistanceDB() { delete distCalculator; }
	
	void generateDB(){}; //doesn't generate a search db 
	void addSequence(Sequence); 
	string getName(int i) { return data[i].getName(); } 
	vector<int> findClosestSequences(Sequence*, int, vector<float>&) const;  // returns indexes of n closest sequences to query
	
private:
	vector<Sequence> data;
	DistCalc* distCalculator;
	
	int templateSeqsLength;
	bool templateAligned;
	
	bool isAligned(string);
	
};

#endif
