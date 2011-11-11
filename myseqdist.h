#ifndef CORRECTDIST_H
#define CORRECTDIST_H


/*
 *  pds.seqdist.h
 *  
 *
 *  Created by Pat Schloss on 8/12/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothurout.h"

/**************************************************************************************************/

class correctDist {
public:
	correctDist(string, int);
	correctDist(int);
	~correctDist(){}
	
	int addSeq(string, string);
	int execute(string);
	
private:
	MothurOut* m;
	int getSequences(string);
	vector<int> fixSequence(string);
	
	int driver(int, int, string);
	int createProcess(string);
	
	double getDist(vector<int>&, vector<int>&);
	int getLastMatch(char, vector<vector<char> >&, int, int, vector<int>&, vector<int>&);
	
	vector<vector<double> > correctMatrix;
	
	vector<vector<int> > sequences;
	
	vector<string> names;	
	int numSeqs;
	int processors;
	vector<int> start;
	vector<int> end;
};

/**************************************************************************************************/

#endif


