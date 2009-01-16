#ifndef SHAREDTHETAYC_H
#define SHAREDTHETAYC_H
/*
 *  sharedthetayc.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedThetaYC estimator on two groups. 
It is a child of the calculator class. */


#include <Carbon/Carbon.h>
#include "calculator.h"

/***********************************************************************/

class SharedThetaYC : public Calculator  {
	
public:
	SharedThetaYC() :  Calculator("SharedThetaYC", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/
 #endif

