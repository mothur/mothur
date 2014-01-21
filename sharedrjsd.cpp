//
//  sharedrjsd.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/21/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "sharedrjsd.h"

/***********************************************************************/
//KLD <- function(x,y) sum(x *log(x/y))
//JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
EstOutput RJSD::getValues(vector<SharedRAbundVector*> shared) {
	try {
        
		data.resize(1,0);
        
        double KLD1 = 0.0;
        double KLD2 = 0.0;
        
        vector<int> countA = shared[0]->getAbundances();
        vector<int> countB = shared[1]->getAbundances();
        double totalA = 0;
        double totalB = 0;
        
		for (int i = 0; i < shared[0]->getNumBins(); i++) {
            totalA += countA[i];
            totalB += countB[i];
        }
        
        for (int i = 0; i < shared[0]->getNumBins(); i++) {
            double tempA = countA[i] / totalA;
            double tempB = countB[i] / totalB;
            
            tempA = countA[i] / totalA;
            tempB = countB[i] / totalB;
            
            if (tempA == 0) { tempA = 0.000001; }
            if (tempB == 0) { tempB = 0.000001; }
            
            double denom = (tempA+tempB)/(double)2.0;
            
            if (tempA != 0) {  KLD1 += tempA * log(tempA/denom); } //KLD(x,m)
            if (tempB != 0) {  KLD2 += tempB * log(tempB/denom); } //KLD(y,m)
            
        }
        
        data[0] = sqrt((0.5*KLD1) + (0.5*KLD2));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "RJSD", "getValues");
		exit(1);
	}
}
