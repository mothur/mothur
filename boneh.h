#ifndef BONEH_H
#define BONEH_H

/*
 *  boneh.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/* This class implements the boneh calculator on single group. 
 It is a child of the calculator class. */

/***********************************************************************/

class Boneh : public Calculator  {
	
public: 
	Boneh(int size) : m(size), Calculator("boneh", 1, false) {};
	EstOutput getValues(SAbundVector*);	
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
private:
	double getV(double, double, double);
	int m;
};


/***********************************************************************/

#endif
