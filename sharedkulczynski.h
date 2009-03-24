#ifndef SHAREDKULCZYNSKI_H
#define SHAREDKULCZYNSKI_H
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

class SharedKulczynski : public Calculator  {
	
public:
	SharedKulczynski() :  Calculator("SharedKulczynski", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/



#endif