#ifndef NSEQS_H
#define NSEQS_H

/*
 *  nseqs.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/16/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class NSeqs : public Calculator {

public:
	NSeqs() : Calculator("nseqs", 1, false) {};
	EstOutput getValues(SAbundVector* rank){
		data.resize(1,0);
		data[0] = (double)rank->getNumSeqs();
		return data;
	}
	EstOutput getValues(vector<SharedRAbundVector*>) {return data;};
};

/***********************************************************************/

#endif
