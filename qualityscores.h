#ifndef QUALITYSCORES
#define QUALITYSCORES

/*
 *  qualityscores.h
 *  Mothur
 *
 *  Created by Pat Schloss on 7/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"
#include "mothurout.h"
#include "sequence.hpp"

/**************************************************************************************************/

class QualityScores {
public:
	QualityScores();
	QualityScores(ifstream&, int);
	string getName();
	void printQScores(ofstream&);
	void trimQScores(int, int);
	void flipQScores();
	bool stripQualThreshold(Sequence&, double);
	bool stripQualRollingAverage(Sequence&, double);
	bool stripQualWindowAverage(Sequence&, int, int, double);
	bool cullQualAverage(Sequence&, double);
private:
	
	double calculateAverage();
	MothurOut* m;
	vector<int> qScores;
	
	string seqName;
	int seqLength;
};
	
/**************************************************************************************************/

#endif
