/*
 *  ssbp.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bergerparker.h"

/***************************************************************/

EstOutput BergerParker::getValues(SAbundVector* rank){
	try {
		data.resize(1,0);
		//Berger-Parker index
		double BP = (double)rank->getMaxRank() / (double)rank->getNumSeqs();
		
		data[0] = BP;
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }

		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "BergerParker", "getValues");
		exit(1);
	}
}

/***********************************************************************/

