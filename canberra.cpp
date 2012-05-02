/*
 *  canberra.cpp
 *  Mothur
 *
 *  Created by westcott on 12/14/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "canberra.h"

/***********************************************************************/

EstOutput Canberra::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		int numSharedOTUS = 0;
		
		double sum = 0.0;
		
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			//is this otu shared
			if ((Aij != 0) && (Bij != 0)) { numSharedOTUS++; }
			
			if ((Aij + Bij) != 0) { 
				sum += ((abs(Aij - Bij)) / (float) (Aij + Bij));
			}
		}
		
		data[0] = (1 / (float) shared[0]->getNumBins()) * sum;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Canberra", "getValues");
		exit(1);
	}
}

/***********************************************************************/
