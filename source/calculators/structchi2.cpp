/*
 *  structchi2.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "structchi2.h"

/***********************************************************************/
EstOutput StructChi2::getValues(vector<SharedRAbundVector*> shared) {
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
				sumOtus[i] += shared[j]->getAbundance(i);
			}
		}
		
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			int A = shared[0]->getAbundance(i);
			int B = shared[1]->getAbundance(i);
			
			double totalTerm = 1 / (float) sumOtus[i];
			double Aterm = A / sumA;
			double Bterm = B / sumB;
			
			sum += (totalTerm * ((Aterm-Bterm)*(Aterm-Bterm)));
		}
				
		data[0] = sqrt((totalSum * sum));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StructChi2", "getValues");
		exit(1);
	}
}
/***********************************************************************/

