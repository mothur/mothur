//
//  shannonrange.cpp
//  Mothur
//
//  Created by SarahsWork on 1/3/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "shannonrange.h"

/***********************************************************************/

EstOutput RangeShannon::getValues(vector<SharedRAbundVector*> shared) {
	try {
        data.resize(3,0);
        
        double commSize = 1e20;
        
        SAbundVector sabund1 = shared[0]->getSAbundVector();
        SAbundVector sabund2 = shared[1]->getSAbundVector();
        
        double sampleSize = 0;
        for (int i = 0; i < sabund1.getNumBins(); i++) {  sampleSize += (sabund1.get(i) * sabund2.get(i));  }
        int aux = ceil(pow((sampleSize+1), 0.33333));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "RangeShannon", "getValues");
		exit(1);
	}
}
/***********************************************************************/
