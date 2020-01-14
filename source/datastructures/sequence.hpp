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

#include "mothurout.h"
#include "utils.hpp"
#include "writer.h"

/**************************************************************************************************/

class Sequence {
    
#ifdef UNIT_TEST
    friend class TestSequence;
#endif

    
public:
    
	Sequence();
	Sequence(string, string);
	Sequence(ifstream&);
    Sequence(ifstream&, string&, bool);
	Sequence(istringstream&);
    #ifdef USE_BOOST
    Sequence(boost::iostreams::filtering_istream&);
    #endif
    ~Sequence() {}
	
	void setName(string);
    string getName();
	void setUnaligned(string);
    string getUnaligned();
	void setAligned(string);
    string getAligned();
    void setComment(string);
    string getComment();
	void setPairwise(string);
	string getPairwise();
	
	string getInlineSeq();
    int getNumNs();
	int getNumBases();
	int getStartPos();
	int getEndPos();
    
    void reverseComplement();
    void trim(int);
	void padToPos(int);
	void padFromPos(int);
    int filterToPos(int); //any character before the pos is changed to . and aligned and unaligned strings changed
    int filterFromPos(int); //any character after the pos is changed to . and aligned and unaligned strings changed
	int getAlignLength();
	int getAmbigBases();
	void removeAmbigBases();
	int getLongHomoPolymer();
    string convert2ints();
    
	void printSequence(ostream&);
    void printSequence(OutputWriter*);
    void printUnAlignedSequence(ostream&);
	
protected:
    
	MothurOut* m;
	void initialize();
	string getSequenceString(ifstream&, int&);
	string getCommentString(ifstream&);
	string getSequenceString(istringstream&, int&);
	string getCommentString(istringstream&);
    string getSequenceName(ifstream&);
    #ifdef USE_BOOST
    string getCommentString(boost::iostreams::filtering_istream&);
    string getSequenceString(boost::iostreams::filtering_istream&, int&);
    string getSequenceName(boost::iostreams::filtering_istream&);
    #endif
    string getSequenceName(istringstream&);
    
    
	string name;
	string unaligned;
	string aligned;
	string pairwise;
    string comment;
	int numBases;
	int alignmentLength;
	int longHomoPolymer;
	int ambigBases;
	int startPos, endPos;
    Utils util;
};

/**************************************************************************************************/

#endif
