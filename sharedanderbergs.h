#ifndef SHAREDANDERBERG_H
#define SHAREDANDERBERG_H
/*
 *  sharedanderberg.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class SharedAnderberg : public Calculator  {
	
	public:
		SharedAnderberg() :  Calculator("SharedAnderberg", 1) {};
		EstOutput getValues(SAbundVector*) {return data;};
		EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
	private:

};

/***********************************************************************/



#endif