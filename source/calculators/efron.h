#ifndef EFRON_H
#define EFRON_H

/*
 *  efron.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/* This class implements the efron calculator on single group. 
 It is a child of the calculator class. */

/***********************************************************************/

class Efron : public Calculator  {
	
public: 
	Efron(int size) : f(size), Calculator("efron", 1, false) {};
	EstOutput getValues(SAbundVector*);	
	EstOutput getValues(vector<RAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Efron"; }
private:
	int f;
};


/***********************************************************************/

#endif
