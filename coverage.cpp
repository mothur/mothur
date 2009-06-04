/*
 *  coverage.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/22/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "coverage.h"

/***********************************************************************/
EstOutput Coverage::getValues(SAbundVector* rank){

	try {
		data.resize(1,0);

		data[0] = 1. - rank->get(1) / (double)rank->getNumSeqs();
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Coverage class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Coverage class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


/***********************************************************************/
