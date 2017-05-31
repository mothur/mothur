#ifndef CANBERRA_H
#define CANBERRA_H

/*
 *  canberra.h
 *  Mothur
 *
 *  Created by westcott on 12/14/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class Canberra : public Calculator  {
	
public:
	Canberra() :  Calculator("canberra", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Canberra"; }
private:
	
};

/***********************************************************************/

#endif

