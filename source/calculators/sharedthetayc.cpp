/*
 *  sharedthetayc.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedthetayc.h"

/***********************************************************************/
EstOutput ThetaYC::getValues(vector<SharedRAbundVector*> shared) {
	try {	
		data.resize(3,0.0000);
		
		double Atotal = (double)shared[0]->getNumSeqs();
		double Btotal = (double)shared[1]->getNumSeqs();
		double thetaYC = 0;
		double pi = 0;
		double qi = 0;
		double a = 0;
		double b = 0;
		double d = 0;
		
		double sumPcubed = 0;
		double sumQcubed = 0;
		double sumPQsq = 0;
		double sumPsqQ = 0;
		
		//calculate the theta denominator sums
		for (int j = 0; j < shared[0]->getNumBins(); j++) {
			//store in temps to avoid multiple repetitive function calls
			pi = shared[0]->get(j) / Atotal;
			qi = shared[1]->get(j) / Btotal;
					
			a += pi * pi;
			b += qi * qi;
			d += pi * qi;
			
			sumPcubed += pi * pi * pi;
			sumQcubed += qi * qi * qi;
			sumPQsq += pi * qi * qi;
			sumPsqQ += pi * pi * qi;
		}

		thetaYC = d / (a + b - d);
		
		if (isnan(thetaYC) || isinf(thetaYC)) { thetaYC = 0; }
		
		double varA = 4 / Atotal * (sumPcubed - a * a);
		double varB = 4 / Btotal * (sumQcubed - b * b);
		double varD = sumPQsq / Atotal + sumPsqQ / Btotal - d * d * (1/Atotal + 1/Btotal);
		double covAD = 2 / Atotal * (sumPsqQ - a * d);
		double covBD = 2 / Btotal * (sumPQsq - b* d);
		
		double varT = d * d * (varA + varB) / pow(a + b - d, (double)4.0) + pow(a+b, (double)2.0) * varD / pow(a+b-d, (double)4.0)
						- 2.0 * (a + b) * d / pow(a + b - d, (double)4.0) * (covAD + covBD);
		
		double ci = 1.95 * sqrt(varT);
		
		data[0] = thetaYC;
		data[1] = thetaYC - ci;
		data[2] = thetaYC + ci;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[2])) { data[2] = 0; }
		
		data[0] = 1.0 - data[0];
        double hold = data[1];
        data[1] = 1.0 - data[2];
        data[2] = 1.0 - hold;
        
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "ThetaYC", "getValues");
		exit(1);
	}
}

/***********************************************************************/
