#ifndef HELLINGER_H
#define HELLINGER_H

/*
 *  hellinger.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class Hellinger : public Calculator  {
	
public:
	Hellinger() :  Calculator("hellinger", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Hellinger"; }
private:
	
};

/***********************************************************************/

#endif


