/*
 *  sharedlennon.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedlennon.h"

/***********************************************************************/

EstOutput Lennon::getValues(vector<SharedRAbundVector*> shared) {
	try {
		double S1, S2, S12, tempA, tempB, min;
		S1 = 0; S2 = 0; S12 = 0; tempA = 0; tempB = 0; min = 0;
		
		/*S1, S2 = number of OTUs observed or estimated in A and B 
		S12=number of OTUs shared between A and B */

		data.resize(1,0);
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared[0]->get(i);
			tempB = shared[1]->get(i);
			
			if (!util.isEqual(tempA, 0)) { S1++; }
			if (!util.isEqual(tempB, 0)) { S2++; }

			//they are shared
			if (!util.isEqual(tempA, 0) && !util.isEqual(tempB, 0)) {	S12++; }
		}
		
		
		tempA = S1 - S12;  tempB = S2 - S12;
		
		if (tempA < tempB) { min = tempA; }
		else { min = tempB; }
		
		data[0] = 1.0 - S12 / (float)(S12 + min);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Lennon", "getValues");
		exit(1);
	}
}

/***********************************************************************/
