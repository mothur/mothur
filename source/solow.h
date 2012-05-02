#ifndef SOLOW_H
#define SOLOW_H

/*
 *  solow.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/* This class implements the solow calculator on single group. 
 It is a child of the calculator class. */

/***********************************************************************/

class Solow : public Calculator  {
	
public: 
	Solow(int size) : f(size), Calculator("solow", 1, false) {};
	EstOutput getValues(SAbundVector*);	
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Solow"; }
private:
	int f;
};


/***********************************************************************/

#endif
