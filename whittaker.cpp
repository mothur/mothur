/*
 *  whittaker.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 4/23/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "whittaker.h"

/***********************************************************************/

EstOutput Whittaker::getValues(vector<SharedRAbundVector*> shared){
	try{
		data.resize(1);

		int countA = 0;
		int countB = 0;
		int sTotal = shared[0]->getNumBins();
		for(int i=0;i<sTotal;i++){
			if(shared[0]->getAbundance(i) != 0){	countA++;	}
			if(shared[1]->getAbundance(i) != 0){	countB++;	}		
		}
		
		data[0] = 2*sTotal/(float)(countA+countB)-1;
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Whittaker class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Whittaker class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
