#ifndef BSTICK_H
#define BSTICK_H
/*
 *  bstick.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the BStick estimator on single group. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class BStick : public Calculator  {
	
public:
	BStick() : Calculator("bstick", 3, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<RAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Bstick"; }
private:
	double invSum(int, double);
	RAbundVector getRAbundVector(SAbundVector*);
	RAbundVector rdata;
};

/***********************************************************************/

#endif
