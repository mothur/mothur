#ifndef KULCZYNSKICODY_H
#define KULCZYNSKICODY_H

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

class KulczynskiCody : public Calculator  {
	
public:
	KulczynskiCody() :  Calculator("KulczynskiCody", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/



#endif