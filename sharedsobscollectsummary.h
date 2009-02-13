#ifndef SHAREDSOBSCOLLECTSUMMARY_H
#define SHAREDSOBSCOLLECTSUMMARY_H

/*
 *  sharedsobscollectsummary.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/12/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This calculator returns the number of shared species between 2 groups. */

#include "calculator.h"

/***********************************************************************/
class SharedSobsCS : public Calculator {

public:
	SharedSobsCS() : Calculator("SharedSobs", 1) {};
	EstOutput getValues(SAbundVector* rank){ return data; };
	EstOutput getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2);
};

/***********************************************************************/

#endif
