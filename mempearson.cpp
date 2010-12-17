/*
 *  mempearson.cpp
 *  Mothur
 *
 *  Created by westcott on 12/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mempearson.h"

/***********************************************************************/
EstOutput MemPearson::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		int nonZeroA = 0;
		int nonZeroB = 0;
		int numOTUS = shared[0]->getNumBins();
		
		//for each otu
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			if (shared[0]->getAbundance(i) != 0) { nonZeroA++; }
			if (shared[1]->getAbundance(i) != 0) { nonZeroB++; }
		}
		
		double numTerm1 = 0.0;
		double numTerm2 = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			int Aij =  shared[0]->getAbundance(i);
			int Bij =  shared[1]->getAbundance(i);
			
			if (Aij > 0) { Aij = 1; }
			if (Bij > 0) { Bij = 1; }
			
			numTerm1 += (Aij - (nonZeroA / (float) numOTUS));
			numTerm2 += (Bij - (nonZeroB / (float) numOTUS));
			
			denomTerm1 += ((Aij - (nonZeroA / (float) numOTUS)) * (Aij - (nonZeroA / (float) numOTUS)));
			denomTerm2 += ((Bij - (nonZeroB / (float) numOTUS)) * (Bij - (nonZeroB / (float) numOTUS)));
		}
		
		denomTerm1 = sqrt(denomTerm1);
		denomTerm2 = sqrt(denomTerm2);
		
		double denom = denomTerm1 * denomTerm2;
		double numerator = numTerm1 * numTerm2;
		
		if (denom != 0) { 
			data[0] = 1.0 - (numerator / denom);
		}else {
			data[0] = 1.0;
		}
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "MemPearson", "getValues");
		exit(1);
	}
}
/***********************************************************************/

