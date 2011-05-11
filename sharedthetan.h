#ifndef THETAN_H
#define THETAN_H
/*
 *  sharedthetan.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedThetaN estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class ThetaN : public Calculator  {
	
public:
	ThetaN() :  Calculator("thetan", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Thetan"; }
private:
	
};

/***********************************************************************/

#endif
