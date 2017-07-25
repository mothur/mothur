/*
 *  spearman.cpp
 *  Mothur
 *
 *  Created by westcott on 12/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "spearman.h"

/***********************************************************************/

EstOutput Spearman::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		
		SAbundVector savA = shared[0]->getSAbundVector();
		SAbundVector savB = shared[1]->getSAbundVector();
		
		double sumRanks = 0.0;
		int numOTUS = shared[0]->getNumBins();
		
		vector<double> rankVectorA(savA.getMaxRank()+1, 0);
		int currentRankA = 0;
		for(int i=savA.getMaxRank();i>0;i--){	
			int numWithAbundanceI = savA.get(i);

			if(numWithAbundanceI > 1){	rankVectorA[i] = (currentRankA + 1 + currentRankA + numWithAbundanceI) / 2.0;	}
			else					{	rankVectorA[i] = currentRankA+numWithAbundanceI;	}
			currentRankA += numWithAbundanceI;
		}
		rankVectorA[0] = (numOTUS + currentRankA + 1) / 2.0;
		
		
		vector<double> rankVectorB(savB.getMaxRank()+1, 0);
		int currentRankB = 0;
		for(int i=savB.getMaxRank();i>0;i--){	
			int numWithAbundanceI = savB.get(i);
			
			if(numWithAbundanceI > 1){	rankVectorB[i] = (currentRankB + 1 + currentRankB + numWithAbundanceI) / 2.0;	}
			else					{	rankVectorB[i] = currentRankB+numWithAbundanceI;	}
			currentRankB += numWithAbundanceI;
		}
		rankVectorB[0] = (numOTUS + currentRankB + 1) / 2.0;
		
		

		for (int i = 0; i < shared[0]->getNumBins(); i++) { 
			int Aij = shared[0]->get(i);
			int Bij = shared[1]->get(i);
			
			float rankA = rankVectorA[Aij];
			float rankB = rankVectorB[Bij];
			
			sumRanks += ((rankA - rankB) * (rankA - rankB));
		}
		data[0] = 1.0 - ((6 * sumRanks) / (float) (numOTUS * ((numOTUS*numOTUS)-1)));
		
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "Spearman", "getValues");
		exit(1);
	}
}

/***********************************************************************/

