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


//Data Structure for a fasta file.

#include "mothur.h"
#include "mothurout.h"

/**************************************************************************************************/

class Sequence {
public:
	Sequence();
	Sequence(string, string);
	Sequence(ifstream&);
	Sequence(istringstream&);
	Sequence(const Sequence& se) : name(se.name), unaligned(se.unaligned), aligned(se.aligned), pairwise(se.pairwise), numBases(se.numBases), startPos(se.startPos), endPos(se.endPos),
									alignmentLength(se.alignmentLength), isAligned(se.isAligned), longHomoPolymer(se.longHomoPolymer), ambigBases(se.ambigBases) { m = MothurOut::getInstance(); }
	
	//these constructors just set the unaligned string to save space
	Sequence(string, string, string);  
	Sequence(ifstream&, string);
	Sequence(istringstream&, string);
	
	void setName(string);
	void setUnaligned(string);
	void setPairwise(string);
	void setAligned(string);
	void setLength();
	void reverseComplement();
	void trim(int);
	
	string convert2ints();
	string getName();
	string getAligned();
	string getPairwise();
	string getUnaligned();
	string getInlineSeq();
	int getNumBases();
	int getStartPos();
	int getEndPos();
	void padToPos(int);
	void padFromPos(int);
	int getAlignLength();
	int getAmbigBases();
	void removeAmbigBases();
	int getLongHomoPolymer();
	bool getIsAligned();
	void printSequence(ostream&);
	
private:
	MothurOut* m;
	void initialize();
	string getSequenceString(ifstream&, int&);
	string getCommentString(ifstream&);
	string getSequenceString(istringstream&, int&);
	string getCommentString(istringstream&);
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
