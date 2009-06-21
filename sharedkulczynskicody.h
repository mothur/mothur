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
	KulczynskiCody() :  Calculator("kulczynskicody", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif
