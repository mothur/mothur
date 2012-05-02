/*
 *  shen.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/18/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "shen.h"
#include "ace.h"

	
/***********************************************************************/	
EstOutput Shen::getValues(SAbundVector* rank){

	try {
		
		data.resize(1,0);
		
		double n = (double)rank->getNumSeqs();
		double f1 = (double)rank->get(1);
		
		Ace* calc = new Ace(abund);
		EstOutput ace = calc->getValues(rank);
		
		double f0 = ace[0]-rank->getNumBins();
		
		data[0] = f0 * (1 - pow(1 - f1/n/f0, f));
		
		delete calc;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Shen", "getValues");
		exit(1);
	}
}


/***********************************************************************/
