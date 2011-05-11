#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H
/*
 *  bootstrap.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the Bootstrap estimator on single group. 
It is a child of the calculator class. */

#include "calculator.h"

/***********************************************************************/


class Bootstrap : public Calculator  {
	
public:
	Bootstrap() : Calculator("bootstrap", 1, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Bootstrap"; }
	
};

/***********************************************************************/

#endif
