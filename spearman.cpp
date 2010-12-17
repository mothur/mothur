/*
 *  spearman.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "spearman.h"

/***********************************************************************/
EstOutput Spearman::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		SAbundVector savA = shared[0]->getSAbundVector();
		SAbundVector savB = shared[1]->getSAbundVector();
		
		double sumRanks = 0.0;
		int numOTUS = shared[0]->getNumBins();
		
		//calc the 2 denominators
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			int rankA = savA.get(Aij);
			int rankB = savB.get(Bij);
			
			sumRanks += ((rankA - rankB) * (rankA - rankB));
		}
		
		data[0] = 1.0 - ((6 * sumRanks) / (float) (numOTUS * (numOTUS-1)));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Soergel", "getValues");
		exit(1);
	}
}
/***********************************************************************/

