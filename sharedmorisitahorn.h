#ifndef SHAREDMORHORN_H
#define SHAREDMORHORN_H
/*
 *  sharedmorisitahorn.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class SharedMorHorn : public Calculator  {
	
public:
	SharedMorHorn() :  Calculator("SharedMorisitaHorn", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/

#endif

