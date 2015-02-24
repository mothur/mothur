/*
 *  manhattan.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "manhattan.h"

/***********************************************************************/
EstOutput Manhattan::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sum = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			sum += abs((Aij - Bij));
		}
		
		data[0] = sum;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Manhattan", "getValues");
		exit(1);
	}
}
/***********************************************************************/
