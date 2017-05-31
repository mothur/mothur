#ifndef MANHATTAN_H
#define MANHATTAN_H

/*
 *  manhattan.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class Manhattan : public Calculator  {
	
public:
	Manhattan() :  Calculator("manhattan", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Manhattan"; }
private:
	
};

/***********************************************************************/

#endif
