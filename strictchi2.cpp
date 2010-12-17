/*
 *  strictchi2.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "strictchi2.h"

/***********************************************************************/
EstOutput StrictChi2::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumA = shared[0]->getNumSeqs();
		double sumB = shared[1]->getNumSeqs();
		double totalSum = 0.0;
		
		for (int i = 0; i < shared.size(); i++) { totalSum += shared[i]->getNumSeqs();  }
		
		vector<int> sumOtus; sumOtus.resize(shared[0]->getNumBins(), 0);
		//for each otu
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			//for each group
			for (int j = 0; j < shared.size(); j++) { 
				
			}
		}
				
		
				
		//data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StrictChi2", "getValues");
		exit(1);
	}
}
/***********************************************************************/
