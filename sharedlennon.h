#ifndef LENNON_H
#define LENNON_H

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

class Lennon : public Calculator  {
	
public:
	Lennon() :  Calculator("Lennon", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/



#endif