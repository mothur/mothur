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
EstOutput SharedBrayCurtis::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
	try {	
		EstOutput data;
		data.resize(1,0);
		
		int sumSharedA, sumSharedB, sumSharedAB, tempA, tempB;
		sumSharedA = 0; sumSharedB = 0; sumSharedAB = 0; 
		
		/*Xi, Yi = abundance of the ith shared OTU in A and B 
		sumSharedA = the sum of all shared otus in A
		sumSharedB = the sum of all shared otus in B
		sumSharedAB = the sum of the minimum otus int all shared otus in AB.
		*/
		
		for (int i = 0; i < shared1->size(); i++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared1->getAbundance(i);
			tempB = shared2->getAbundance(i);

			
			if ((tempA != 0) && (tempB != 0)) {//they are shared
				sumSharedA += tempA;
				sumSharedB += tempB;
				
				//sum the min of tempA and tempB
				if (tempA < tempB) { sumSharedAB += tempA; }
				else  { sumSharedAB += tempB; }				
			}
		}
		
		data[0] = (2 * sumSharedAB) / (float)( sumSharedA + sumSharedB);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
				
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedBrayCurtis class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedBrayCurtis class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	


}

/***********************************************************************/
