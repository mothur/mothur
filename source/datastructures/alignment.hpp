#ifndef DPALIGNMENT_H
#define DPALIGNMENT_H

/*
 *  dpalignment.h
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *  This is a class for an abstract datatype for classes that implement various types of alignment	algorithms.
 *	As of 12/18/08 these included alignments based on blastn, needleman-wunsch, and the	Gotoh algorithms
 * 
 */

#include "mothur.h"
#include "alignmentcell.hpp"
#include "currentfile.h"

/**************************************************************************************************/

class Alignment {
	
public:
	Alignment(int);
    Alignment(int, int);
	Alignment();
	virtual ~Alignment();
	virtual void align(string, string) = 0;
    virtual void alignPrimer(string, string) {}
	
	
//	float getAlignmentScore();
	string getSeqAAln();
	string getSeqBAln();
    map<int, int> getSeqAAlnBaseMap();
    map<int, int> getSeqBAlnBaseMap();
	int getCandidateStartPos();
	int getCandidateEndPos();
	int getTemplateStartPos();
	int getTemplateEndPos();
	
	int getPairwiseLength();
	void resize(int);
	int getnRows() { return nRows; }
//	int getLongestTemplateGap();

protected:
	void traceBack();
	string seqA, seqAaln;
	string seqB, seqBaln;
	int seqAstart, seqAend;
	int seqBstart, seqBend;
	int pairwiseLength;
	int nRows, nCols, lA, lB;
	vector<vector<AlignmentCell> > alignment;
    map<int, int> ABaseMap;
    map<int, int> BBaseMap;
	MothurOut* m;
    CurrentFile* current;
};

/**************************************************************************************************/

#endif
