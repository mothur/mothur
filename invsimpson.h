#ifndef INVSIMPSON
#define INVSIMPSON

/*
 *  invsimpson.h
 *  Mothur
 *
 *  Created by Pat Schloss on 8/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class InvSimpson : public Calculator  {
	
public:
	InvSimpson() : Calculator("invsimpson", 3, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
};

/***********************************************************************/

#endif
