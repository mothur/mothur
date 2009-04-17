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

/* This class implements the LogSD estimator on single group. 
It is a child of the calculator class. */

/***********************************************************************/

class Geom : public Calculator  {
	
public:
	Geom() : Calculator("geom", 3) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};

private:
	double kEq(double, double);
	RAbundVector getRAbundVector(SAbundVector*);
	RAbundVector rdata;
};

/***********************************************************************/

#endif





