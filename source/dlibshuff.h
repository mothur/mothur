#ifndef DLIBSHUFF
#define DLIBSHUFF

/*
 *  dlibshuff.h
 *  Mothur
 *
 *  Created by Pat Schloss on 4/8/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "fullmatrix.h"
#include "libshuff.h"

class DLibshuff : public Libshuff {

public:
	DLibshuff(FullMatrix*, int, float, float);
	vector<vector<double> > evaluateAll();
	float evaluatePair(int, int);

private:
	int numDXs;
	double dCalculate(int, int);
	vector<int> calcN(vector<double>);
};

#endif
