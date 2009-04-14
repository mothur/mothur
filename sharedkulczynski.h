#ifndef KULCZYNSKI_H
#define KULCZYNSKI_H
/*
 *  sharedkulczynski.h
 *  Mothur
 *
 *  Created by John Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class Kulczynski : public Calculator  {
	
public:
	Kulczynski() :  Calculator("Kulczynski", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/



#endif