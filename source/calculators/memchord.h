#ifndef MEMCHORD_H
#define MEMCHORD_H

/*
 *  memchord.h
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/***********************************************************************/

class MemChord : public Calculator  {
	
public:
	MemChord() :  Calculator("memchord", 1, false) {};  
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Memchord"; }
private:
	
};

/***********************************************************************/

#endif

