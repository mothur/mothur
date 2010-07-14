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
		
		float Atotal = 0;
		float Btotal = 0;
		float thetaYC = 0;
		float pi = 0;
		float qi = 0;
		float a = 0;
		float b = 0;
		float d = 0;
		
		float sumPcubed = 0;
		float sumQcubed = 0;
		float sumPQsq = 0;
		float sumPsqQ = 0;
		
		//get the total values we need to calculate the theta denominator sums
		for (int i = 0; i < shared[0]->size(); i++) {
			//store in temps to avoid multiple repetitive function calls
			Atotal += (float)shared[0]->getAbundance(i);
			Btotal += (float)shared[1]->getAbundance(i);
		}
		
		//calculate the theta denominator sums
		for (int j = 0; j < shared[0]->size(); j++) {
			//store in temps to avoid multiple repetitive function calls
			pi = shared[0]->getAbundance(j) / Atotal;
			qi = shared[1]->getAbundance(j) / Btotal;
					
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
		
		float varA = 4 / Atotal * (sumPcubed - a * a);
		float varB = 4 / Btotal * (sumQcubed - b * b);
		float varD = sumPQsq / Atotal + sumPsqQ / Btotal - d * d * (1/Atotal + 1/Btotal);
		float covAD = 2 / Atotal * (sumPsqQ - a * d);
		float covBD = 2 / Btotal * (sumPQsq - b* d);
		
		float varT = d * d * (varA + varB) / pow(a + b - d, (float)4.0) + pow(a+b, (float)2.0) * varD / pow(a+b-d, (float)4.0)
						- 2.0 * (a + b) * d / pow(a + b - d, (float)4.0) * (covAD + covBD);
		
		float ci = 1.95 * sqrt(varT);
		
		data[0] = thetaYC;
		data[1] = thetaYC - ci;
		data[2] = thetaYC + ci;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "ThetaYC", "getValues");
		exit(1);
	}
}

/***********************************************************************/
