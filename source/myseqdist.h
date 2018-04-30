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
#include "utils.hpp"

/**************************************************************************************************/

class correctDist {
public:
	correctDist(string, int);
	correctDist(int);
	~correctDist(){}
	
	int addSeq(string, string);
	void execute(string);
	
private:
	MothurOut* m;
    Utils util;
    vector<vector<int> > sequences;
    vector<string> names;
    int processors;
    
	int getSequences(string);
	vector<int> fixSequence(string);
	int createProcess(string);
};

/**************************************************************************************************/

#endif


