#ifndef THETAN_H
#define THETAN_H
/*
 *  sharedthetan.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedThetaN estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class ThetaN : public Calculator  {
	
public:
	ThetaN() :  Calculator("ThetaN", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/

#endif
