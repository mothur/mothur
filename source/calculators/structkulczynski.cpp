/*
 *  structkulczynski.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "structkulczynski.h"

/***********************************************************************/
EstOutput StructKulczynski::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumA = 0.0;
		double sumB = 0.0;
		double sumMin = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int A = shared[0]->getAbundance(i);
			int B = shared[1]->getAbundance(i);
			
			sumA += A;
			sumB += B;
			sumMin += min(A, B);
		}
		
		data[0] = 1.0 - (0.5 * ((sumMin / sumA) + (sumMin / sumB)));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StructKulczynski", "getValues");
		exit(1);
	}
}
/***********************************************************************/


