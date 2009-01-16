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

#include <string>


/**************************************************************************************************/

class Sequence {
public:
	Sequence();
	Sequence(ifstream&);
	void setName(string);
	void setUnaligned(string);
	void setPairwise(string);
	void setAligned(string);
	string convert2ints();
	string getSeqName();
	string getAligned();
	string getPairwise();
	string getUnaligned();
private:
	string name;
	string unaligned;
	string pairwise;
	string aligned;
};

/**************************************************************************************************/

#endif
