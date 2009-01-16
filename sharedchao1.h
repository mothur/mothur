#ifndef SHAREDCHAO1_H
#define SHAREDCHAO1_H
/*
 *  sharedchao1.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the Sharedchao1 estimator on two groups. 
It is a child of the calculator class. */


#include <Carbon/Carbon.h>
#include "calculator.h"

/***********************************************************************/

class SharedChao1 : public Calculator  {
	
public: 
	SharedChao1() : Calculator("SharedChao", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
};

/***********************************************************************/
#endif