/*
 *  hellinger.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "hellinger.h"

/***********************************************************************/
EstOutput Hellinger::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumA = 0.0;
		double sumB = 0.0;
		
		//calc the 2 denominators
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			sumA += shared[0]->getAbundance(i);
			sumB += shared[1]->getAbundance(i);
		}
		
		
		//calc sum
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			double term1 = sqrt((Aij / sumA));
			double term2 = sqrt((Bij / sumB));
			
			sum += ((term1 - term2) * (term1 - term2));
		}
		
		data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Hellinger", "getValues");
		exit(1);
	}
}
/***********************************************************************/


