#ifndef SHAREDSOBS_H
#define SHAREDSOBS_H
/*
 *  sharedsobs.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedSobs estimator on two groups for the shared rarefaction command. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/
class SharedSobs : public Calculator {

public:
	SharedSobs() : Calculator("sharedsobs", 1, false) {};
	EstOutput getValues(SAbundVector* rank){ return data; };
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/sharedsobs"; }
};

/***********************************************************************/

#endif

