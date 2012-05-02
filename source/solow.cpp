/*
 *  solow.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "solow.h"
#include <math.h>

	
/***********************************************************************/	
EstOutput Solow::getValues(SAbundVector* rank){

	try {
		data.resize(1,0);
		
		double n = (double)rank->getNumSeqs();
		double f1 = (double)rank->get(1);
		double f2 = (double)rank->get(2);

		data[0] = f1*f1/2/f2 * (1 - pow(1 - 2*f2/n/f1, f));
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Solow", "getValues");
		exit(1);
	}
}


/***********************************************************************/
