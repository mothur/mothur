#ifndef GOODSCOVERAGE_H
#define GOODSCOVERAGE_H

/*
 *  goodscoverage.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/* This class implements the LogSD estimator on single group. 
It is a child of the calculator class. */

/***********************************************************************/

class GoodsCoverage : public Calculator  {
	
public:
	GoodsCoverage() : Calculator("goodscoverage", 1, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<RAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/GoodsCoverage"; }

private:
};

/***********************************************************************/

#endif
