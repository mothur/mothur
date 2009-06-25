/*
 *  sharedsorest.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedsorest.h"
#include "chao1.h"
#include "sharedchao1.h"

/***********************************************************************/

EstOutput SorEst::getValues(vector<SharedRAbundVector*> shared) {
	try {
		EstOutput S1, S2, S12;
		S12.resize(1,0);
		S1.resize(3,0);
		S2.resize(3,0); 
		
		/*S1, S2 = number of OTUs estimated in A and B using the Chao estimator
		S12 = estimated number of OTUs shared between A and B using the SharedChao estimator*/

		data.resize(1,0);
		
		SharedChao1* sharedChao = new SharedChao1();
		Chao1* chaoS1 = new Chao1();
		Chao1* chaoS2 = new Chao1();
		SAbundVector* chaoS1Sabund = new SAbundVector();
		SAbundVector* chaoS2Sabund = new SAbundVector();
		
		*chaoS1Sabund = shared[0]->getSAbundVector();
		*chaoS2Sabund = shared[1]->getSAbundVector();
		
		S12 = sharedChao->getValues(shared);
		S1 = chaoS1->getValues(chaoS1Sabund);
		S2 = chaoS2->getValues(chaoS2Sabund);
		
		data[0] = (2 * S12[0]) / (float)(S1[0] + S2[0]);
		 
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		delete sharedChao;
		delete chaoS1;
		delete chaoS2;
		delete chaoS1Sabund;
		delete chaoS2Sabund;
		
		return data;
	}
	catch(exception& e) {
		errorOut(e, "SorEst", "getValues");
		exit(1);
	}
}

/***********************************************************************/
