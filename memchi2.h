#ifndef MEMCHI2_H
#define MEMCHI2_H

/*
 *  memchi2.h
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class MemChi2 : public Calculator  {
	
public:
	MemChi2() :  Calculator("memchi2", 1, false, true) {};  //the true means this calculator needs all groups to calculate the pair value
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Memchi2"; }
private:
	
};

/***********************************************************************/

#endif
