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
//cout << "singles = " << singles << " doubles = " << doubles << " sobs = " << sobs << endl;	
		double chao = sobs + singles*(singles-1)/(2*(doubles+1));
	
		if(singles > 0 && doubles > 0){
			chaovar = singles*(singles-1)/(2*(doubles+1))
					+ singles*pow(2*singles-1, 2)/(4*pow(doubles+1,2))
					+ pow(singles, 2)*doubles*pow(singles-1, 2)/(4*pow(doubles+1,4));
		}
		else if(singles > 0 && doubles == 0){
			chaovar = singles*(singles-1)/2 + singles*pow(2*singles-1, 2)/4 - pow(singles, 4)/(4*chao); 
		}
		else if(singles == 0){
			chaovar = sobs*exp(-1*rank->getNumSeqs()/sobs)*(1-exp(-1*rank->getNumSeqs()/sobs));
		}
			
		double chaohci, chaolci;
	
		if(singles>0){
			double denom = pow(chao-sobs,2);
			double c = exp(1.96*pow((log(1+chaovar/denom)),0.5));
			chaolci = sobs+(chao-sobs)/c;//chao lci
			chaohci = sobs+(chao-sobs)*c;//chao hci
		}
		else{
			double p = exp(-1*rank->getNumSeqs()/sobs);
			chaolci = sobs/(1-p)-1.96*pow(sobs*p/(1-p), 0.5);
			chaohci = sobs/(1-p)+1.96*pow(sobs*p/(1-p), 0.5);
			if(chaolci < sobs){	chaolci = sobs;	}
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
		errorOut(e, "Chao1", "getValues");
		exit(1);
	}
}

/***********************************************************************/
