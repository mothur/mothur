/*
 *  efron.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "efron.h"

/***********************************************************************/
EstOutput Efron::getValues(SAbundVector* rank){

	try {
		data.resize(1,0);
		
		double n = (double)rank->getNumSeqs();
		if(f > n || f == 0) {	f = n;	}
		
		double sum = 0;
		for(int i = 1; i < rank->size(); i++){
			sum += pow(-1., i+1) * pow(((double)f / n), i) * (double)(rank->get(i));
		}
		data[0] = sum;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Efron", "getValues");
		exit(1);
	}
}


/***********************************************************************/
