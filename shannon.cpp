/*
 *  shannon.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "shannon.h"

/***********************************************************************/

EstOutput Shannon::getValues(SAbundVector* rank){
	try {
		//vector<double> shannonData(3,0);
		data.resize(3,0);
	
		double shannon = 0.0000;  //hprime
		double hvara=0.0000;
    
		double maxRank = rank->getMaxRank();
		int sampled = rank->getNumSeqs();
		int sobs = rank->getNumBins();
	
		for(int i=1;i<=maxRank;i++){
			double p = ((double) i)/((double)sampled);
			shannon += (double)rank->get(i)*p*log(p); //hprime
			hvara  += (double)rank->get(i)*p*pow(log(p),2);
		}
		shannon = -shannon;
    
		double hvar = (hvara-pow(shannon,2))/(double)sampled+(double)sobs/(double)(2*sampled*sampled);
    
		double ci = 0;
	
		if(hvar>0){
			ci = 1.96*pow(hvar,0.5);
		}
	
		double shannonhci = shannon + ci;
		double shannonlci = shannon - ci;
		
		data[0] = shannon;
		data[1] = shannonlci;
		data[2] = shannonhci;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }

		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Shannon class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Shannon class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
