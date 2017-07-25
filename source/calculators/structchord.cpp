/*
 *  structchord.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "structchord.h"

/***********************************************************************/
EstOutput StructChord::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		double sumAj2 = 0.0;
		double sumBj2 = 0.0;
		
		//calc the 2 denominators
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->get(i);
			int Bij = shared[1]->get(i);
			
			//(Aij) ^ 2
			sumAj2 += (Aij * Aij);
			sumBj2 += (Bij * Bij);
		}
		
		sumAj2 = sqrt(sumAj2);
		sumBj2 = sqrt(sumBj2);
		
		//calc sum
		double sum = 0.0;
		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			
			int Aij = shared[0]->get(i);
			int Bij = shared[1]->get(i);
			
			sum += (((Aij / sumAj2) - (Bij / sumBj2)) * ((Aij / sumAj2) - (Bij / sumBj2)));
		}
		
		data[0] = sqrt(sum);
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "StructChord", "getValues");
		exit(1);
	}
}
/***********************************************************************/
