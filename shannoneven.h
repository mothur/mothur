#ifndef SHANNONEVEN
#define SHANNONEVEN

/*
 *  shannoneven.h
 *  Mothur
 *
 *  Created by Pat Schloss on 8/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class ShannonEven : public Calculator  {
	
public:
	ShannonEven() : Calculator("shannoneven", 1, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
};

/***********************************************************************/

#endif
