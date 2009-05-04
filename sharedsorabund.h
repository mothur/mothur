#ifndef SORABUND_H
#define SORABUND_H
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


#include "calculator.h"

/***********************************************************************/

class SorAbund : public Calculator  {
	
public:
	SorAbund() :  Calculator("SorAbund", 3, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	UVEst* uv;
	
};
/***********************************************************************/

#endif
