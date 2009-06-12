#ifndef SHEN_H
#define SHEN_H

/*
 *  shen.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/18/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"
#include "globaldata.hpp"

/* This class implements the shen calculator on single group. 
 It is a child of the calculator class. */

/***********************************************************************/

class Shen : public Calculator  {
	
public: 
	Shen(int size, int n) : m(size), abund(n), Calculator("shen", 1, false) {};
	EstOutput getValues(SAbundVector*);	
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
private:
	int m;
	int abund;
	GlobalData* globaldata;
};


/***********************************************************************/

#endif
