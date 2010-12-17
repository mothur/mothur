#ifndef STRICTCHORD_H
#define STRICTCHORD_H

/*
 *  strictchord.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class StrictChord : public Calculator  {
	
public:
	StrictChord() :  Calculator("strictchord", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
private:
	
};

/***********************************************************************/

#endif



