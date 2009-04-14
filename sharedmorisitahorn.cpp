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
EstOutput MorHorn::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
	try {	
		data.resize(1,0);
		
		int Atotal, Btotal, tempA, tempB;
		Atotal = 0; Btotal = 0; 
		float  morhorn, sumSharedA, sumSharedB, a, b, d;
		morhorn = 0.0; sumSharedA = 0.0; sumSharedB = 0.0; a = 0.0; b = 0.0; d = 0.0;
		
		//get the total values we need to calculate the theta denominator sums
		for (int i = 0; i < shared1->size(); i++) {
			//store in temps to avoid multiple repetitive function calls
			Atotal += shared1->getAbundance(i);
			Btotal += shared2->getAbundance(i);
		}
		
		//calculate the theta denominator sums
		for (int j = 0; j < shared1->size(); j++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared1->getAbundance(j);
			tempB = shared2->getAbundance(j);
			
			//they are shared
			if ((tempA != 0) && (tempB != 0)) {
				if (Atotal != 0)	{ sumSharedA = (tempA / (float)Atotal); }
				if (Btotal != 0)	{ sumSharedB = (tempB / (float)Btotal); }
			
				a += sumSharedA * sumSharedA;
				b += sumSharedB * sumSharedB;
				d += sumSharedA * sumSharedB;
			}
		}

		morhorn = (2 * d) / (float) (a + b);
		
		if (isnan(morhorn) || isinf(morhorn)) { morhorn = 0; }
		
		data[0] = morhorn;
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the MorHorn class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the MorHorn class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/