#ifndef COVERAGE_H
#define COVERAGE_H

/*
 *  coverage.h
 *  Mothur
 *
 *  Created by Pat Schloss on 4/22/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "calculator.h"

/* This class implements the coverage estimator on single group. 
 It is a child of the calculator class. */

/***********************************************************************/

class Coverage : public Calculator  {
	
public: 
	Coverage() : Calculator("coverage", 1, false) {};
	EstOutput getValues(SAbundVector*);	
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Coverage"; }
};


/***********************************************************************/

#endif
