#ifndef STRICTEUCLIDEAN_H
#define STRICTEUCLIDEAN_H

/*
 *  stricteuclidean.h
 *  Mothur
 *
 *  Created by westcott on 12/14/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class StrictEuclidean : public Calculator  {
	
public:
	StrictEuclidean() :  Calculator("stricteuclidean", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif


