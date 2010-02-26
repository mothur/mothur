/*
 *  simpson.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "simpson.h"

/***********************************************************************/

EstOutput Simpson::getValues(SAbundVector* rank){
	try {
		//vector<double> simpsonData(3,0);
		data.resize(3,0);
		double simpson = 0.0000;
		double ci = 0;
	
		double maxRank = (double)rank->getMaxRank();
		double sampled = (double)rank->getNumSeqs();
		double sobs = (double)rank->getNumBins();
	
		double firstTerm = 0;
		double secondTerm = 0;
	
		if(sobs != 0){
			double simnum=0.0000;
		
			for(int i=1;i<=maxRank;i++){
				simnum += (double)(rank->get(i)*i*(i-1));
			}
			
			simpson = simnum / (sampled*(sampled-1));
		
			for(int i=1;i<=maxRank;i++){
				double piI = (double) i / (double)sampled;
				firstTerm += rank->get(i) * pow(piI, 3);
				secondTerm += rank->get(i) * pow(piI, 2);
			}
		
			double var = (4.0 / sampled) * (firstTerm - secondTerm*secondTerm);
			ci = 1.95 * pow(var, 0.5);
		}
	
		double simpsonlci = simpson - ci;
		double simpsonhci = simpson + ci;
		
		data[0] = simpson;
		data[1] = simpsonlci;
		data[2] = simpsonhci;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
	
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Simpson", "getValues");
		exit(1);
	}
}

/***********************************************************************/
