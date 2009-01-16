#ifndef SHAREDSORABUND_H
#define SHAREDSORABUND_H
/*
 *  sharedsorabund.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedSorAbund estimator on two groups. 
It is a child of the calculator class. */


#include <Carbon/Carbon.h>
#include "calculator.h"

/***********************************************************************/

class SharedSorAbund : public Calculator  {
	
public:
	SharedSorAbund() :  Calculator("SharedSorAbund", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	UVEst* uv;
	
};
/***********************************************************************/

#endif
