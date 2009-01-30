#ifndef SHAREDJABUND_H
#define SHAREDJABUND_H
/*
 *  sharedjabund.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedJAbund estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class SharedJAbund : public Calculator  {
	
public:
	SharedJAbund() :  Calculator("SharedJAbund", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	UVEst* uv;
	
};

/***********************************************************************/

#endif
