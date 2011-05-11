#ifndef NPSHANNON_H
#define NPSHANNON_H

/*
 *  npshannon.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the NPShannon estimator on single group. 
It is a child of the calculator class. */

#include "calculator.h"

/***********************************************************************/

class NPShannon : public Calculator  {
	
public:
	NPShannon() : Calculator("npshannon", 1, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Npshannon"; }
private:
	
};

/***********************************************************************/

#endif
