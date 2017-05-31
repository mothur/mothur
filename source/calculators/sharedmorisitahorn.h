#ifndef MORHORN_H
#define MORHORN_H
/*
 *  sharedmorisitahorn.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class MorHorn : public Calculator  {
	
public:
	MorHorn() :  Calculator("morisitahorn", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Morisitahorn"; }
private:
	
};

/***********************************************************************/

#endif

