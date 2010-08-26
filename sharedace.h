#ifndef SHAREDACE_H
#define SHAREDACE_H
/*
 *  sharedace.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedAce estimator on two groups. 
It is a child of the calculator class. */

#include "calculator.h"

/***********************************************************************/

class SharedAce : public Calculator  {
	
public:
	SharedAce(int n=10) : abund(n),  Calculator("sharedace", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	int abund;
};

/***********************************************************************/

#endif
