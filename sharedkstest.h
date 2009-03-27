#ifndef SHAREDKSTEST_H
#define SHAREDKSTEST_H
/*
 *  kstest.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the KSTest estimator on 2 groups. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class SharedKSTest : public Calculator  {
	
public:
	SharedKSTest() : Calculator("sharedkstest", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
};

/***********************************************************************/

#endif
