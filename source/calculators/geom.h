#ifndef GEOM_H
#define GEOM_H

/*
 *  geom.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/* This class implements the geometric estimator on single group. 
It is a child of the calculator class. */

/***********************************************************************/

class Geom : public Calculator  {
	
public:

	Geom() : Calculator("geometric", 3, false) {};

	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Geometric"; }
private:
	double kEq(double, double);
	RAbundVector getRAbundVector(SAbundVector*);
	RAbundVector rdata;
};

/***********************************************************************/

#endif





