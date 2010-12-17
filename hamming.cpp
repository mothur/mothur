/*
 *  hamming.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "hamming.h"

/***********************************************************************/
EstOutput Hamming::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		int numA = 0;
		int numB = 0;
		int numShared = 0;
		
		//calc the 2 denominators
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			int A = shared[0]->getAbundance(i);
			int B = shared[1]->getAbundance(i);
			
			if (A != 0) { numA++; }
			if (B != 0) { numB++; }
			if ((A != 0) && (B != 0)) { numShared++; }
		}
		
		data[0] = numA + numB - (2 * numShared);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Hamming", "getValues");
		exit(1);
	}
}
/***********************************************************************/

