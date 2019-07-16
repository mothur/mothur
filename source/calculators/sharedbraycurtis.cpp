/*
 *  sharedbraycurtis.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedbraycurtis.h"

/***********************************************************************/
//This is used by SharedJAbund and SharedSorAbund
EstOutput BrayCurtis::getValues(vector<SharedRAbundVector*> shared) {
	try {	
		data.resize(1,0);
		
		double sumSharedA, sumSharedB, sumSharedAB, tempA, tempB;
		sumSharedA = shared[0]->getNumSeqs(); sumSharedB = shared[1]->getNumSeqs(); sumSharedAB = 0;
		
		/*Xi, Yi = abundance of the ith shared OTU in A and B 
		sumSharedA = the number of otus in A
		sumSharedB = the sum of all shared otus in B
		sumSharedAB = the sum of the minimum otus int all shared otus in AB.
		*/
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared[0]->get(i);
			tempB = shared[1]->get(i);
				
			//sum the min of tempA and tempB
			if (tempA < tempB) { sumSharedAB += tempA; }
			else  { sumSharedAB += tempB; }				
		}
		
		data[0] = 1.0 - (2 * sumSharedAB) / (float)( sumSharedA + sumSharedB);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
				
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "BrayCurtis", "getValues");
		exit(1);
	}
}

/***********************************************************************/
