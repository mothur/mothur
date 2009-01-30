#ifndef SHAREDJCLASS_H
#define SHAREDJCLASS_H
/*
 *  sharedjclass.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedJclass estimator on two groups. 
It is a child of the calculator class. */

#include "calculator.h"

/***********************************************************************/

class SharedJclass : public Calculator  {
	
public:
	SharedJclass() :  Calculator("SharedJclass", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/

#endif
