#ifndef SEQUENCE_H
#define SEQUENCE_H

/*
 *  sequence.h
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	A sequence object has three components: i) an accession number / name, ii) the unaligned primary sequence, iii) a
 *	pairwise aligned sequence, and iv) a sequence that is aligned to a reference alignment.  This class has methods
 *	to set and get these values for the other classes where they are needed.  *
 *
 */
using namespace std;

#include "mothur.h"


/**************************************************************************************************/

class Sequence {
public:
	Sequence();
	Sequence(string, string);
	Sequence(ifstream&);
	
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
	int getLength();  //the greater of the lengths of unaligned and aligned
	int getUnalignLength();
	int getAlignLength();
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
