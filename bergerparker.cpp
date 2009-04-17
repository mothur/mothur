/*
 *  ssbp.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bergerparker.h"

/***************************************************************/

EstOutput BergerParker::getValues(SAbundVector* rank){
	try {
		data.resize(1,0);
		//Berger-Parker index
		double BP = (double)rank->getNumSeqs()/(double)rank->getMaxRank();
		//cout << "BP index = " << 1/BP << "\n\n";
		
		data[0] = BP;
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }

		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the NPShannon class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the NPShannon class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

