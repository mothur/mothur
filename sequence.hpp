#ifndef SEQUENCE_H
#define SEQUENCE_H

/*
 *  sequence.h
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */
using namespace std;

#include "mothur.h"


/**************************************************************************************************/

class Sequence {
public:
	Sequence();
	Sequence(string, string);
	
	void setName(string);
	void setUnaligned(string);
	void setPairwise(string);
	void setAligned(string);
	void setLength();
	
	string convert2ints();
	string getName();
	string getAligned();
	string getPairwise();
	string getUnaligned();
	int getLength();
	void printSequence(ostream&);
	
private:
	string name;
	string unaligned;
	string aligned;
	string pairwise;
	int length;
	int lengthAligned;
};

/**************************************************************************************************/

#endif
