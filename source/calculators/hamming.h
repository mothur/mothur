#ifndef HAMMING_H
#define HAMMING_H

/*
 *  hamming.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class Hamming : public Calculator  {
	
public:
	Hamming() :  Calculator("hamming", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Hamming"; }
private:
	
};

/***********************************************************************/

#endif

