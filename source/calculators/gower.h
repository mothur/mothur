#ifndef GOWER_H
#define GOWER_H

/*
 *  gower.h
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class Gower : public Calculator  {
	
public:
	Gower() :  Calculator("gower", 1, false, true) {};  //the true means this calculator needs all groups to calculate the pair value
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Gower"; }
private:
	
};

/***********************************************************************/

#endif

