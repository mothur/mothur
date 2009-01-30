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
	Ace(int n=10) : abund(n), Calculator("ACE", 3) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};
private:
	int abund;
};

/***********************************************************************/

#endif
