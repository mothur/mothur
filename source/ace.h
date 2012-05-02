#ifndef ACE_H
#define ACE_H

/*
 *  ace.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the Ace estimator on a single group. 
It is a child of the calculator class. */

#include "calculator.h"

/***********************************************************************/

class Ace : public Calculator  {
	
public:
	Ace(int n) : abund(n), Calculator("ace", 3, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Ace"; }
private:
	int abund;
};

/***********************************************************************/

#endif
