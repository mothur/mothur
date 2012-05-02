/*
 *  smithwilson.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 8/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "smithwilson.h"

/***********************************************************************/

EstOutput SmithWilson::getValues(SAbundVector* rank){
	try {

		data.resize(1,0);

		double maxRank = rank->getMaxRank();
		double sobs = rank->getNumBins();
		
		double innerSum = 0;
		for(int i=1;i<=maxRank;i++){
			innerSum += rank->get(i) * log(i);
		}
		innerSum /= sobs;

		double outerSum = 0;
		for(int i=1;i<=maxRank;i++){
			outerSum += rank->get(i) * (log(i) - innerSum) * (log(i) - innerSum);
		}
		outerSum /= sobs;

		if(outerSum > 0){
			data[0] = 1.0000 - 2.0000 / (3.14159 * atan(outerSum));
		}
		else{
			data[0] = 1.0000;
		}
		
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "InvSimpson", "getValues");
		exit(1);
	}
}

/***********************************************************************/
