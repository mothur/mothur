/*
 *  sharedsobscollectsummary.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/12/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedsobscollectsummary.h"

/***********************************************************************/
//This returns the number of shared species observed in several groups.  
//The shared vector is each groups sharedrabundvector.

EstOutput SharedSobsCS::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	try {
		data.resize(1,0);
		int observed = 0;
		int tempA, tempB;

		//loop through the species in each group
		for (int k = 0; k < shared1->size(); k++) {
			tempA = shared1->getAbundance(k); //store in temps to avoid calling getAbundance multiple times
			tempB = shared2->getAbundance(k);

			//if you have found a new species
			if ((tempA != 0) && (tempB != 0)) {//they are shared
				observed++;
			}
		}

		data[0] = observed;
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSobsCS class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSobsCS class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/