/*
 *  speciesprofile.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "speciesprofile.h"

/***********************************************************************/
EstOutput SpeciesProfile::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumA = 0.0;
		double sumB = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			sumA += shared[0]->get(i);
			sumB += shared[1]->get(i);
		}
		
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			int A = shared[0]->get(i);
			int B = shared[1]->get(i);
			
			sum += (((A / sumA) - (B / sumB)) * ((A / sumA) - (B / sumB)));
		}
		
		data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "SpeciesProfile", "getValues");
		exit(1);
	}
}
/***********************************************************************/

