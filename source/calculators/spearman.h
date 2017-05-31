#ifndef SPEARMAN_H
#define SPEARMAN_H

/*
 *  spearman.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class Spearman : public Calculator  {
	
public:
	Spearman() :  Calculator("spearman", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Spearman"; }
private:
	
};

/***********************************************************************/

#endif





