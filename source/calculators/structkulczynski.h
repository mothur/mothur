#ifndef STRUCTKULCZYNSKI_H
#define STRUCTKULCZYNSKI_H

/*
 *  structkulczynski.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class StructKulczynski : public Calculator  {
	
public:
	StructKulczynski() :  Calculator("structkulczynski", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Structkulczynski"; }
private:
	
};

/***********************************************************************/

#endif



