#ifndef SHAREDMARCZEWSKI_H
#define SHAREDMARCZEWSKI_H

/*
 *  sharedmarczewski.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/8/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class SharedMarczewski : public Calculator  {
	
public:
	SharedMarczewski() : Calculator("sharedmarczewski", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
};

/***********************************************************************/

#endif