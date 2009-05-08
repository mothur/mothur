#ifndef SIMPSON_H
#define SIMPSON_H
/*
 *  simpson.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the Simpson estimator on single group. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class Simpson : public Calculator  {

public:
	Simpson() : Calculator("Simpson", 3, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
};

/***********************************************************************/

#endif
