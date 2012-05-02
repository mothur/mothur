/*
 *  memeuclidean.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "memeuclidean.h"

/***********************************************************************/
EstOutput MemEuclidean::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			int A = shared[0]->getAbundance(i);
			int B = shared[1]->getAbundance(i);
			
			if (A > 0) { A = 1; }
			if (B > 0) { B = 1; }
			
			sum += ((A-B)*(A-B));
		}
		
		data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "MemEuclidean", "getValues");
		exit(1);
	}
}
/***********************************************************************/

