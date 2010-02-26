/*
 *  coverage.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/22/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "coverage.h"

/***********************************************************************/
EstOutput Coverage::getValues(SAbundVector* rank){

	try {
		data.resize(1,0);

		data[0] = 1. - rank->get(1) / (double)rank->getNumSeqs();
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Coverage", "getValues");
		exit(1);
	}
}


/***********************************************************************/
