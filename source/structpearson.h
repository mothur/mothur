#ifndef STRUCTPEARSON_H
#define STRUCTPEARSON_H

/*
 *  structpearson.h
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class StructPearson : public Calculator  {
	
public:
	StructPearson() :  Calculator("structpearson", 1, false) {};
	EstOutput getValues(SAbundVector*) {return data;};
	EstOutput getValues(vector<SharedRAbundVector*>);
	string getCitation() { return "http://www.mothur.org/wiki/Structpearson"; }
private:
	
};

/***********************************************************************/

#endif


