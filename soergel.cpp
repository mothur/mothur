/*
 *  soergel.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "soergel.h"

/***********************************************************************/
EstOutput Soergel::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumNum = 0.0;
		double sumMax = 0.0;
		
		//calc the 2 denominators
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->getAbundance(i);
			int Bij = shared[1]->getAbundance(i);
			
			sumNum += abs((Aij - Bij));
			sumMax += max(Aij, Bij);
		}
		
		data[0] = sumNum / sumMax;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Soergel", "getValues");
		exit(1);
	}
}
/***********************************************************************/

