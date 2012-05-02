/*
 *  sharedsobs.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedsobs.h"

/***********************************************************************/
//This returns the number of unique species observed in several groups.  
//The shared vector is each groups sharedrabundvector.

EstOutput SharedSobs::getValues(vector<SharedRAbundVector*> shared){
	try {
		data.resize(1,0);
		double observed = 0;

		//loop through the species in each group
		for (int k = 0; k < shared[0]->getNumBins(); k++) {
			//if you have found a new species
			if (shared[0]->getAbundance(k) != 0) { observed++; } 
			else if ((shared[0]->getAbundance(k) == 0) && (shared[1]->getAbundance(k) != 0)) { observed++; }
		}

		data[0] = observed;
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedSobs", "getValues");
		exit(1);
	}
}

/***********************************************************************/
