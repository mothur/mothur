#ifndef STRUCTCHORD_H
#define STRUCTCHORD_H

/*
 *  structchord.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class StructChord : public Calculator  {
	
public:
	StructChord() :  Calculator("structchord", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Structchord"; }
private:
	
};

/***********************************************************************/

#endif



