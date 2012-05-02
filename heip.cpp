/*
 *  heip.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 8/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "heip.h"
#include "shannon.h"

/***********************************************************************/

EstOutput Heip::getValues(SAbundVector* rank){
	try {
		data.resize(1,0.0000);
		vector<double> shanData(3,0);
		Shannon* shannon = new Shannon();
		shanData = shannon->getValues(rank);
		
		long int sobs = rank->getNumBins();
		if(sobs > 1){
			data[0] = (exp(shanData[0])-1) / (sobs - 1);;
		}
		else{
			data[0] = 1;
		}
		
		delete shannon;
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Heip", "getValues");
		exit(1);
	}
}

/***********************************************************************/
