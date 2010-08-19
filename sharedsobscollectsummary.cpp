/*
 *  sharedsobscollectsummary.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/12/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedsobscollectsummary.h"

/***********************************************************************/
//This returns the number of shared species observed in several groups.  
//The shared vector is each groups sharedrabundvector.

EstOutput SharedSobsCS::getValues(vector<SharedRAbundVector*> shared){
	try {
		data.resize(1,0);
		double observed = 0;
		int numGroups = shared.size();

		for (int i = 0; i < shared[0]->size(); i++) {
			//get bin values and set sharedByAll 
			bool sharedByAll = true;
			for (int j = 0; j < numGroups; j++) {
				if (shared[j]->getAbundance(i) == 0) { sharedByAll = false; }
			}
			
			//they are shared
			if (sharedByAll == true) {  observed++;  }
		}

		data[0] = observed;
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedSobsCS", "getValues");
		exit(1);
	}
}

/***********************************************************************/
