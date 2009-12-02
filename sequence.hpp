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
	void reverseComplement();
	
	string convert2ints();
	string getName();
	string getAligned();
	string getPairwise();
	string getUnaligned();
	int getNumBases();
	int getStartPos();
	int getEndPos();
	int getAlignLength();
	int getAmbigBases();
	int getLongHomoPolymer();
	bool getIsAligned();
	void printSequence(ostream&);
	
private:
	void initialize();
	string getSequenceString(ifstream&);
	string getCommentString(ifstream&);
	string name;
	string unaligned;
	string aligned;
	string pairwise;
	int numBases;
	int alignmentLength;
	bool isAligned;
	int longHomoPolymer;
	int ambigBases;
	int startPos, endPos;
};

/**************************************************************************************************/

#endif
