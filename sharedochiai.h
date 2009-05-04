#ifndef OCHIAI_H
#define OCHIAI_H
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

class Ochiai : public Calculator  {
	
public:
	Ochiai() :  Calculator("Ochiai", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/


#endif