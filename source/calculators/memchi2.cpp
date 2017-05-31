/*
 *  memchi2.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "memchi2.h"

/***********************************************************************/
EstOutput MemChi2::getValues(vector<RAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		int nonZeroA = 0;
		int nonZeroB = 0;
		int totalOtus = shared[0]->getNumBins();
		//int totalGroups = shared.size();
		
		//for each otu
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			if (shared[0]->get(i) != 0) { nonZeroA++; }
			if (shared[1]->get(i) != 0) { nonZeroB++; }
		}
		
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			int A = shared[0]->get(i);
			int B = shared[1]->get(i);
			
			if (A > 0) { A = 1; }
			if (B > 0) { B = 1; }
			
			double Aterm = A / (float) nonZeroA;
			double Bterm = B / (float) nonZeroB;

			int incidence = 0;
			for(int j=0;j<shared.size();j++){
				if(shared[j]->get(i) != 0){	incidence++;	}
			}
			
			if(incidence != 0){
				sum += (((Aterm-Bterm)*(Aterm-Bterm))/incidence);
			}
		}
		
		data[0] = sqrt(totalOtus * sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "MemChi2", "getValues");
		exit(1);
	}
}
/***********************************************************************/

