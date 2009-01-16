#ifndef SIMPSON_H
#define SIMPSON_H
/*
 *  simpson.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab Umass Amherst. All rights reserved.
 *
 */

/* This class implements the Simpson estimator on single group. 
It is a child of the calculator class. */


#include <Carbon/Carbon.h>
#include "calculator.h"

/***********************************************************************/

class Simpson : public Calculator  {

public:
	Simpson() : Calculator("Simpson", 3) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};
};

/***********************************************************************/

#endif