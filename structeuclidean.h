#ifndef STRUCTEUCLIDEAN_H
#define STRUCTEUCLIDEAN_H

/*
 *  structeuclidean.h
 *  Mothur
 *
 *  Created by westcott on 12/14/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class StructEuclidean : public Calculator  {
	
public:
	StructEuclidean() :  Calculator("structeuclidean", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Structeuclidean"; }
private:
	
};

/***********************************************************************/

#endif


