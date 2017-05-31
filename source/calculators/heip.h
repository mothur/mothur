#ifndef HEIP
#define HEIP

/*
 *  heip.h
 *  Mothur
 *
 *  Created by Pat Schloss on 8/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class Heip : public Calculator  {
	
public:
	Heip() : Calculator("heip", 1, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<RAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Heip"; }
};

/***********************************************************************/

#endif
