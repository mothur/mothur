/*
 *  structeuclidean.cpp
 *  Mothur
 *
 *  Created by westcott on 12/14/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "structeuclidean.h"

/***********************************************************************/
EstOutput StructEuclidean::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sum = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			//(Aij - Bij) ^ 2
			sum += ((Aij - Bij) * (Aij - Bij));
		}
		
		data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StructEuclidean", "getValues");
		exit(1);
	}
}
/***********************************************************************/

