//
//  sharedjsd.cpp
//  Mothur
//
//  Created by SarahsWork on 12/9/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "sharedjsd.h"

/***********************************************************************/
//KLD <- function(x,y) sum(x *log(x/y))
//JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
EstOutput JSD::getValues(vector<SharedRAbundVector*> shared) {
	try {
        
		data.resize(1,0);
        
        double KLD1 = 0.0;
        double KLD2 = 0.0;

        double totalA = shared[0]->getNumSeqs();
        double totalB = shared[1]->getNumSeqs();
        
        for (int i = 0; i < shared[0]->getNumBins(); i++) {
            double tempA = shared[0]->get(i) / totalA;
            double tempB = shared[1]->get(i) / totalB;
            
            if (util.isEqual(tempA, 0)) { tempA = 0.000001; }
            if (util.isEqual(tempB, 0)) { tempB = 0.000001; }

            double denom = (tempA+tempB)/(double)2.0;

            if (!util.isEqual(tempA, 0)) {  KLD1 += tempA * log(tempA/denom); } //KLD(x,m)
            if (!util.isEqual(tempB, 0)) {  KLD2 += tempB * log(tempB/denom); } //KLD(y,m)
        }

            
        data[0] = ((0.5*KLD1) + (0.5*KLD2));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "JSD", "getValues");
		exit(1);
	}
}

/***********************************************************************/
