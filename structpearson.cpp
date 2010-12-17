/*
 *  structpearson.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "structpearson.h"

/***********************************************************************/
EstOutput StructPearson::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumA = 0.0;
		double sumB = 0.0;
		int numOTUS = shared[0]->getNumBins();
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			sumA += shared[0]->getAbundance(i);
			sumB += shared[1]->getAbundance(i);
		}
		
		double numTerm1 = 0.0;
		double numTerm2 = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			int Aij =  shared[0]->getAbundance(i);
			int Bij =  shared[1]->getAbundance(i);
			
			numTerm1 += (Aij - (sumA / (float) numOTUS));
			numTerm2 += (Bij - (sumB / (float) numOTUS));
			
			denomTerm1 += ((Aij - (sumA / (float) numOTUS)) * (Aij - (sumA / (float) numOTUS)));
			denomTerm2 += ((Bij - (sumB / (float) numOTUS)) * (Bij - (sumB / (float) numOTUS)));
		}
		
		denomTerm1 = sqrt(denomTerm1);
		denomTerm2 = sqrt(denomTerm2);
		
		double denom = denomTerm1 * denomTerm2;
		double numerator = numTerm1 * numTerm2;
		
		data[0] = 1.0 - (numerator / denom);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StructPearson", "getValues");
		exit(1);
	}
}
/***********************************************************************/
