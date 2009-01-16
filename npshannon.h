#ifndef NPSHANNON_H
#define NPSHANNON_H

/*
 *  npshannon.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the NPShannon estimator on single group. 
It is a child of the calculator class. */

#include <Carbon/Carbon.h>
#include "calculator.h"

/***********************************************************************/

class NPShannon : public Calculator  {
	
public:
	NPShannon() : Calculator("NPShannon", 1) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};
private:
	
};

/***********************************************************************/

#endif