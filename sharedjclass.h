#ifndef JCLASS_H
#define JCLASS_H
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

class Jclass : public Calculator  {
	
public:
	Jclass() :  Calculator("jclass", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif
