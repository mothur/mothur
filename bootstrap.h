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

#include <Carbon/Carbon.h>
#include "calculator.h"

/***********************************************************************/


class Bootstrap : public Calculator  {
	
public:
	Bootstrap() : Calculator("Bootstrap", 1) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};
	
};

/***********************************************************************/

#endif