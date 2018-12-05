/*
 *  gower.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "gower.h"

/***********************************************************************/
EstOutput Gower::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		vector<int> maxOtus; maxOtus.resize(shared[0]->getNumBins());
		vector<int> minOtus; minOtus.resize(shared[0]->getNumBins());
		//for each otu
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			
			//set otus min and max to first one
			minOtus[i] = shared[0]->get(i);
			maxOtus[i] = shared[0]->get(i);
			
			//for each group
			for (int j = 1; j < shared.size(); j++) { 
				maxOtus[i] = max(shared[j]->get(i), maxOtus[i]);
				minOtus[i] = min(int(shared[j]->get(i)), minOtus[i]);
			}
		}
		
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			int A = shared[0]->get(i);
			int B = shared[1]->get(i);
			
			double numerator = abs(A - B);
			double denominator = maxOtus[i] - minOtus[i];
				
			if (!util.isEqual(denominator, 0)) { sum += (numerator / denominator); }
		}
		
		data[0] = sum;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Gower", "getValues");
		exit(1);
	}
}
/***********************************************************************/

