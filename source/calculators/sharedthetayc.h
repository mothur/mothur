#ifndef THETAYC_H
#define THETAYC_H
/*
 *  sharedthetayc.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class implements the SharedThetaYC estimator on two groups. 
It is a child of the calculator class. */


#include "calculator.h"

/***********************************************************************/

class ThetaYC : public Calculator  {
	
public:
	ThetaYC() :  Calculator("thetayc", 3, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Thetayc"; }
private:
	
};

/***********************************************************************/

#endif

