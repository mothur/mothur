/*
 *  bootstrap.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bootstrap.h"

/***********************************************************************/

EstOutput Bootstrap::getValues(SAbundVector* rank){
	try {
		//vector<double> bootData(3,0);
		data.resize(1,0);
		double maxRank = (double)rank->getMaxRank();
		double sampled = rank->getNumSeqs();
		double sobs = rank->getNumBins();

		double boot = (double)sobs;

		for(int i=1;i<=maxRank;i++){
			boot += (double)rank->get(i)*pow((1.0-(double)i/(double)sampled),sampled);
		}
		
		data[0] = boot;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
	
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Bootstrap", "getValues");
		exit(1);
	}
}

/***********************************************************************/
