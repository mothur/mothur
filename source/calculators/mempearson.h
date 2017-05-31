#ifndef MEMPEARSON_H
#define MEMPEARSON_H

/*
 *  mempearson.h
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class MemPearson : public Calculator  {
	
public:
	MemPearson() :  Calculator("mempearson", 1, false) {};  
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<RAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Mempearson"; }
private:
	
};

/***********************************************************************/

#endif
