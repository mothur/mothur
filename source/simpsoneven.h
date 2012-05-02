#ifndef SIMPSONEVEN_H
#define SIMPSONEVEN_H

/*
 *  simpsoneven.h
 *  Mothur
 *
 *  Created by Pat Schloss on 8/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/

class SimpsonEven : public Calculator  {
	
public:
	SimpsonEven() : Calculator("simpsoneven", 1, false) {};
	EstOutput getValues(SAbundVector*);
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
	string getCitation() { return "http://www.mothur.org/wiki/Simpsoneven"; }
};

/***********************************************************************/

#endif

