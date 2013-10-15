/*
 *  sharedjest.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "sharedchao1.h"
#include "chao1.h"
#include "sharedjest.h"


/***********************************************************************/

EstOutput Jest::getValues(vector<SharedRAbundVector*> shared) {
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
        
        //chaoS1Sabund->print(cout);
        //chaoS2Sabund->print(cout);
		
		S12 = sharedChao->getValues(shared);
		S1 = chaoS1->getValues(chaoS1Sabund);
		S2 = chaoS2->getValues(chaoS2Sabund);
        
        //cout << S12[0] << '\t' << S1[0] << '\t' << S2[0] << endl;
		
		data[0] = 1.0 - S12[0] / (float)(S1[0] + S2[0] - S12[0]);
        //cout << data[0] << endl;
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		delete sharedChao;
		delete chaoS1;
		delete chaoS2;
		delete chaoS1Sabund;
		delete chaoS2Sabund;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Jest", "getValues");
		exit(1);
	}
}

/***********************************************************************/
