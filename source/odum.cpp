/*
 *  odum.cpp
 *  Mothur
 *
 *  Created by westcott on 12/14/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "odum.h"

/***********************************************************************/

EstOutput Odum::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumNum = 0.0;
		double sumDenom = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			sumNum += abs(Aij - Bij);
			sumDenom += (Aij + Bij);
		}
		
		data[0] = sumNum / sumDenom;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Odum", "getValues");
		exit(1);
	}
}

/***********************************************************************/
