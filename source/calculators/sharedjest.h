#ifndef JEST_H
#define JEST_H
/*
 *  sharedjest.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedJest estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class Jest : public Calculator  {
	
public:
	Jest() :  Calculator("jest", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Jest"; }
private:
	
};

/***********************************************************************/

#endif
