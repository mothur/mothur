#ifndef SHAREDMARCZEWSKI_H
#define SHAREDMARCZEWSKI_H

/*
 *  sharedmarczewski.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class SharedMarczewski : public Calculator  {
	
public:
	SharedMarczewski() : Calculator("sharedmarczewski", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Sharedmarczewski"; }
private:
};

/***********************************************************************/

#endif
