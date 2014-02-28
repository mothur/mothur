//
//  shannonrange.cpp
//  Mothur
//
//  Created by SarahsWork on 1/3/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "shannonrange.h"

/***********************************************************************/

EstOutput RangeShannon::getValues(SAbundVector* rank){
	try {
        data.resize(3,0);
        
        double commSize = 1e20;
        double sampleSize = rank->getNumSeqs();
        
        double aux = ceil(pow(sampleSize+1, (1/(double)3)));
        double est0 = rank->get(1)+1;
        if (aux > est0) { est0 = aux; } //est0 = max(rank->get(1)+1, aux)
        
        vector<double> ests;
        double numr = 0.0;
        for (int i = 1; i < rank->getNumBins()-1; i++) {
            
            if (m->control_pressed) { break; }
            
            int abund = rank->get(i);
            
            if (abund != 0) {
            
                int abundNext = rank->get(i+1);
                if (abundNext == 0) {  numr = aux;  }
                else {
                    if (abundNext+1 > aux) { numr = abundNext+1; } //numr = max(abundNext+1, aux)
                    else { numr = aux; }
                }
                double denr = aux;
                if (abund > aux) { denr = abund; } //denr = max(abund, aux)
                ests.push_back((abund+1)*numr/denr);
             }
        }
        numr = aux;
        
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "RangeShannon", "getValues");
		exit(1);
	}
}
/***********************************************************************/
