/*
 *  simpsoneven.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 8/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "simpsoneven.h"
#include "invsimpson.h"

/***********************************************************************/

EstOutput SimpsonEven::getValues(SAbundVector* rank){
	try {
		data.resize(1,0);

		InvSimpson* simp = new InvSimpson();
		vector<double> invSimpData = simp->getValues(rank);
		
		data[0] = invSimpData[0] / double(rank->getNumBins());
		
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "SimpsonEven", "getValues");
		exit(1);
	}
}

/***********************************************************************/

