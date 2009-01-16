/*
 *  chao1.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chao1.h"

/***********************************************************************/
EstOutput Chao1::getValues(SAbundVector* rank){
	try {
		data.resize(3,0);
	
		double sobs = (double)rank->getNumBins();
	
		double singles = (double)rank->get(1);
		double doubles = (double)rank->get(2);
		double chaovar = 0.0000;
	
		double chao = sobs + pow(singles,2)/(2*(doubles+1)) - (singles*doubles/(2*pow(doubles+1,2)));
	
		if(doubles==0){chaovar=0;}
		else{
			double g=singles/doubles;
			chaovar = doubles*(0.25*pow(g,4)+pow(g,3)+0.5*pow(g,2));
		}

	
		double chaohci, chaolci;
	
		if(chao==sobs){
			double ci = 1.96*pow(chaovar,0.5);
			chaolci = chao-ci;//chao lci
			chaohci = chao+ci;//chao hci
		}
		else{
			double denom = pow(chao-sobs,2);
			double c = exp(1.96*pow((log(1+chaovar/denom)),0.5));
			chaolci = sobs+(chao-sobs)/c;//chao lci
			chaohci = sobs+(chao-sobs)*c;//chao hci
		}
		
		data[0] = chao;
		data[1] = chaolci;
		data[2] = chaohci;

	    if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Chao1 class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Chao1 class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
};

/***********************************************************************/
