#ifndef BRAYCURTIS_H
#define BRAYCURTIS_H
/*
 *  sharedbraycurtis.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include "calculator.h"

/***********************************************************************/

class BrayCurtis : public Calculator  {
	
public:
	BrayCurtis() :  Calculator("BrayCurtis", 1) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(SharedRAbundVector*, SharedRAbundVector*);
private:
	
};

/***********************************************************************/



#endif
