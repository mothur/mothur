#ifndef BSTICK_H
#define BSTICK_H
/*
 *  bstick.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the BStick estimator on single group. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class BStick : public Calculator  {
	
public:
	BStick() : Calculator("bstick", 3) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};

private:
	double invSum(int, double);
	RAbundVector getRAbundVector(SAbundVector*);
	RAbundVector rdata;
};

/***********************************************************************/

#endif
