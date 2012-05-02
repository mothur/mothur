/*
 *  sharedsorclass.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedsorclass.h"

/***********************************************************************/

EstOutput SorClass::getValues(vector<SharedRAbundVector*> shared) {
	try {
		double S1, S2, S12, tempA, tempB;
		S1 = 0; S2 = 0; S12 = 0; tempA = 0; tempB = 0; 
		
		/*S1, S2 = number of OTUs observed or estimated in A and B 
		S12=number of OTUs shared between A and B */

		data.resize(1,0);
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared[0]->getAbundance(i);
			tempB = shared[1]->getAbundance(i);

			//find number of bins in shared1 and shared2
			if (tempA != 0) { S1++; }
			if (tempB != 0) { S2++; } 
			
			//they are shared
			if ((tempA != 0) && (tempB != 0)) {	S12++; }
		}
		
		data[0] = 1.0-(2 * S12) / (float)(S1 + S2);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "SorClass", "getValues");
		exit(1);
	}
}

/***********************************************************************/
