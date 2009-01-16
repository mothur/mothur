/*
 *  sharedsobs.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedsobs.h"

/***********************************************************************/
//This returns the number of unique species observed in several groups.  
//The shared vector is each groups sharedrabundvector.

EstOutput SharedSobs::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	try {
		data.resize(1,0);
		int observed = 0;;

		//loop through the species in each group
		for (int k = 0; k < shared1->size(); k++) {
			//if you have found a new species
			if (shared1->getAbundance(k) != 0) { observed++; } 
			else if ((shared1->getAbundance(k) == 0) && (shared2->getAbundance(k) != 0)) { observed++; }
		}

		data[0] = observed;
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSobs class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSobs class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
