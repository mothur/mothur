#ifndef BERGERPARKER_H
#define BERGERPARKER_H
/*
 *  bergerparker.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/*This class implements the SSBP estimator on single group. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class BergerParker : public Calculator  {
	
public:
	BergerParker() : Calculator("bergerparker", 1) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*) {return data;};

private:
};

/***********************************************************************/

#endif

