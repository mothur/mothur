#ifndef SHAREDOCHIAI_H
#define SHAREDOCHIAI_H
/*
 *  sharedochiai.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class SharedOchiai : public Calculator  {
	
public:
	SharedOchiai() :  Calculator("SharedOchiai", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/


#endif