#ifndef SHANNON_H
#define SHANNON_H
/*
 *  shannon.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the Shannon estimator on single group. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class Shannon : public Calculator  {
	
public:
	Shannon() : Calculator("shannon", 3, false) {};
	EstOutput getValues(SAbundVector* rank);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
private:
	
};

/***********************************************************************/

#endif
