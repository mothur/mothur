#ifndef SOREST_H
#define SOREST_H
/*
 *  sharedsorest.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedSorEst estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class SorEst : public Calculator  {
	
public:
	SorEst() :  Calculator("sorest", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Sorest"; }
private:
	
};

/***********************************************************************/

#endif
