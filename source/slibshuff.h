#ifndef SLIBSHUFF
#define SLIBSHUFF

/*
 *  slibshuff.h
 *  Mothur
 *
 *  Created by Pat Schloss on 4/8/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "fullmatrix.h"
#include "libshuff.h"

class SLibshuff : public Libshuff {

public:
	SLibshuff(FullMatrix*, int, float);
	vector<vector<double> > evaluateAll();
	float evaluatePair(int, int);
	
private:
	double sCalculate(int, int);
};

#endif
