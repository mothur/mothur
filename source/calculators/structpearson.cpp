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
		
		int numOTUS = shared[0]->getNumBins();
		double averageA = shared[0]->getNumSeqs() /  (float) numOTUS;
		double averageB = shared[1]->getNumSeqs() / (float) numOTUS;
		
		double numTerm = 0.0;
		double denomTerm1 = 0.0;
		double denomTerm2 = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			int Aij =  shared[0]->get(i);
			int Bij =  shared[1]->get(i);
			
			
			numTerm += ((Aij - averageA) * (Bij - averageB));
			denomTerm1 += ((Aij - averageA) * (Aij - averageA));
			denomTerm2 += ((Bij - averageB) * (Bij - averageB));
		}
		
		denomTerm1 = sqrt(denomTerm1);
		denomTerm2 = sqrt(denomTerm2);
		
		double denom = denomTerm1 * denomTerm2;
		
		data[0] = (numTerm / denom);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StructPearson", "getValues");
		exit(1);
	}
}
/***********************************************************************/
