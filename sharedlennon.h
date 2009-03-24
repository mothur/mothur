#ifndef SHAREDLENNON_H
#define SHAREDLENNON_H

/*
 *  sharedlennon.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class SharedLennon : public Calculator  {
	
public:
	SharedLennon() :  Calculator("SharedLennon", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/



#endif