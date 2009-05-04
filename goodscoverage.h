#ifndef GOODSCOVERAGE_H
#define GOODSCOVERAGE_H

/*
 *  goodscoverage.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/8/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
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
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};

private:
};

/***********************************************************************/

#endif
