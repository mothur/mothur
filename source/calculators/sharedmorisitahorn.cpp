/*
 *  sharedmorisitahorn.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedmorisitahorn.h"

/***********************************************************************/
EstOutput MorHorn::getValues(vector<SharedRAbundVector*> shared) {
	try {	
		data.resize(1,0);
		
		double Atotal, Btotal, tempA, tempB;
		Atotal = 0; Btotal = 0; 
		double  morhorn, sumSharedA, sumSharedB, a, b, d;
		morhorn = 0.0; sumSharedA = 0.0; sumSharedB = 0.0; a = 0.0; b = 0.0; d = 0.0;
		
		//get the total values we need to calculate the theta denominator sums
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			//store in temps to avoid multiple repetitive function calls
			Atotal += shared[0]->get(i);
			Btotal += shared[1]->get(i);
		}
		
		//calculate the denominator sums
		for (int j = 0; j < shared[0]->getNumBins(); j++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared[0]->get(j);
			tempB = shared[1]->get(j);
			float relA = tempA / Atotal;
			float relB = tempB / Btotal;
			
			a += relA * relA;
			b += relB * relB;
			d += relA * relB;
		}

		morhorn = 1- (2 * d) / (a + b);

		if (isnan(morhorn) || isinf(morhorn)) { morhorn = 1; }
		
		data[0] = morhorn;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "MorHorn", "getValues");
		exit(1);
	}
}

/***********************************************************************/
