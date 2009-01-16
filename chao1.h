#ifndef CHAO1_H
#define CHAO1_H
/*
 *  chao1.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include "calculator.h"

/* This class implements the Ace estimator on single group. 
It is a child of the calculator class. */

/***********************************************************************/

class Chao1 : public Calculator  {
	
public: 
	Chao1() : Calculator("Chao", 3) {};
	EstOutput getValues(SAbundVector*);	
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};
};


/***********************************************************************/

#endif