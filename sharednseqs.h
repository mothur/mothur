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
	SharedNSeqs() : Calculator("SharedNSeqs", 1) {};
	EstOutput getValues(SAbundVector* rank){ return data; };
	EstOutput getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
		data.resize(1,0);
		data[0] = (double)shared1->getNumSeqs() + (double)shared2->getNumSeqs();
		return data;
	}
};

/***********************************************************************/

#endif
