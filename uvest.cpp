/*
 *  uvest.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "uvest.h"

/***********************************************************************/
//This is used by SharedJAbund and SharedSorAbund
EstOutput UVEst::getUVest(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
	try {	
		EstOutput results;
		results.resize(2,0);
		
		int S12, Atotal, Btotal, f1A, f2A, f1B, f2B, sumSharedA, sumSharedB, sumSharedA1, sumSharedB1, tempA, tempB;
		S12 = 0; Atotal = 0; Btotal = 0; f1A = 0; f2A = 0; f1B = 0; f2B = 0; sumSharedA = 0; sumSharedB = 0; sumSharedA1 = 0; sumSharedB1 = 0;
		float Upart1, Upart2, Upart3, Vpart1, Vpart2, Vpart3, Uest, Vest;
		Upart1 = 0.0; Upart2 = 0.0; Upart3 = 0.0; Vpart1 = 0.0; Vpart2 = 0.0; Vpart3 = 0.0;
		
		/*Xi, Yi = abundance of the ith shared OTU in A and B 
		ntotal, mtotal = total number of sequences sampled in A and B 
		I(•) = if the argument, •, is true then I(•) is 1; otherwise it is 0. 
		sumSharedA = the sum of all shared otus in A
		sumSharedB = the sum of all shared otus in B
		sumSharedA1 = the sum of all shared otus in A where B = 1
		sumSharedB1 = the sum of all shared otus in B where A = 1 */
		
		for (int i = 0; i < shared1->size(); i++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared1->getAbundance(i);
			tempB = shared2->getAbundance(i);

			Atotal += tempA;
			Btotal += tempB;
			
			if ((tempA != 0) && (tempB != 0)) {//they are shared
				sumSharedA += tempA;
				sumSharedB += tempB;
				
				//does A have one or two
				if (tempA == 1)			{ f1A++; sumSharedB1 += tempB;}
				else if (tempA == 2)	{ f2A++; }
				
				//does B have one or two
				if (tempB == 1)			{ f1B++; sumSharedA1 += tempA;}
				else if (tempB == 2)	{ f2B++; }
			}
		}
		
		Upart1 = sumSharedA / (float) Atotal;
		Upart2 = ((Btotal - 1) * f1B) / (float) (Btotal * 2 * f2B);
		Upart3 = sumSharedA1 / (float) Atotal;
		
		if (isnan(Upart1) || isinf(Upart1)) { Upart1 = 0; }
		if (isnan(Upart2) || isinf(Upart2)) { Upart2 = 0; }
		if (isnan(Upart3) || isinf(Upart3)) { Upart3 = 0; }
		
		Uest = Upart1 + (Upart2 * Upart3);
		
		Vpart1 = sumSharedB / (float) Btotal;
		Vpart2 = ((Atotal - 1) * f1A) / (float) (Atotal * 2 * f2A);
		Vpart3 = sumSharedB1 / (float) Btotal;
		
		if (isnan(Vpart1) || isinf(Vpart1)) { Vpart1 = 0; }
		if (isnan(Vpart2) || isinf(Vpart2)) { Vpart2 = 0; }
		if (isnan(Vpart3) || isinf(Vpart3)) { Vpart3 = 0; }
		
		Vest = Vpart1 + (Vpart2 * Vpart3);
		
		if (Uest > 1) { Uest = 1; }
		if (Vest > 1) { Vest = 1; }
		
		results[0] = Uest;
		results[1] = Vest;
		
		return results;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UVEst class Function getUVest. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UVEst class Function getUVest. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	


}

/***********************************************************************/
