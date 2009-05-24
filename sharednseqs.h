#ifndef SHAREDNSEQS_H
#define SHAREDNSEQS_H

/*
 *  sharednseqs.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/16/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/***********************************************************************/
class SharedNSeqs : public Calculator {

public:
	SharedNSeqs() : Calculator("sharednseqs", 1, false) {};
	EstOutput getValues(SAbundVector* rank){ return data; };
	EstOutput getValues(vector<SharedRAbundVector*> shared) {
		data.resize(1,0);
		data[0] = (double)shared[0]->getNumSeqs() + (double)shared[1]->getNumSeqs();
		return data;
	}
};

/***********************************************************************/

#endif
