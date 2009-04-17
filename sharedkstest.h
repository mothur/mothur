#ifndef KSTEST_H
#define KSTEST_H
/*
 *  kstest.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/*This class implements the KSTest estimator on 2 groups. 
It is a child of the calculator class.*/ 

/***********************************************************************/

class KSTest : public Calculator  {
	
public:
	KSTest() : Calculator("kstest", 3) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
};

/***********************************************************************/

#endif
