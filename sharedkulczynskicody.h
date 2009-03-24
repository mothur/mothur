#ifndef SHAREDKULCZYNSKICODY_H
#define SHAREDKULCZYNSKICODY_H

/*
 *  sharedkulczynskicody.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class SharedKulczynskiCody : public Calculator  {
	
public:
	SharedKulczynskiCody() :  Calculator("SharedKulczynskiCody", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/



#endif