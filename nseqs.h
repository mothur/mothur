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
	NSeqs() : Calculator("NSeqs", 1) {};
	EstOutput getValues(SAbundVector* rank){
		data.resize(1,0);
		data[0] = (double)rank->getNumSeqs();
		return data;
	}
	EstOutput getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {return data;};
};

/***********************************************************************/

#endif
