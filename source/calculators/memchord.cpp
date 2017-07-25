/*
 *  memchord.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "memchord.h"

/***********************************************************************/
EstOutput MemChord::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double nonZeroA = 0;
		double nonZeroB = 0;
		
		//for each otu
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			if (shared[0]->get(i) != 0) { nonZeroA++; }
			if (shared[1]->get(i) != 0) { nonZeroB++; }
		}
		
		nonZeroA = sqrt(nonZeroA);
		nonZeroB = sqrt(nonZeroB);
		
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			int A = shared[0]->get(i);
			int B = shared[1]->get(i);
			
			if (A > 0) { A = 1; }
			if (B > 0) { B = 1; }
			
			double Aterm = A / nonZeroA;
			double Bterm = B / nonZeroB;
			
			sum += ((Aterm-Bterm)*(Aterm-Bterm));
		}
		
		data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "MemChord", "getValues");
		exit(1);
	}
}
/***********************************************************************/

