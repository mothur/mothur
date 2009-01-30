#ifndef RARECALC_H
#define RARECALC_H
/*
 *  rarecalc.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class is not currently used by Mothur */

#include "calculator.h"

/***********************************************************************/
class RareCalc {

public:
	RareCalc(RAbundVector* b) : bins(b), numSeqs(b->getNumSeqs()), maxRank(b->getMaxRank()), numBins(b->getNumBins()) {	bMatrix = binomial(numSeqs+1);	};
	EstOutput getValues(int);
	string getName()	{	return "RareCalc";	}
private:
	RAbundVector* bins;
	vector<vector<double> > bMatrix;
	int numSeqs, maxRank, numBins;
};

/***********************************************************************/

#endif
