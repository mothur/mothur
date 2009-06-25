
/*
 *  sharedace.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedace.h"

/***********************************************************************/

EstOutput SharedAce::getValues(vector<SharedRAbundVector*> shared) {
	try {
		data.resize(1,0);
		string label;
		label = shared[0]->getLabel();

		int fARare, fBRare, S12Rare, S12Abund, S12, f11, tempA, tempB, t10, t01, t11, t21, t12, t22, C12Numerator;
		fARare = 0; fBRare = 0; S12Rare = 0; S12Abund = 0; S12 = 0; f11 = 0; t10 = 0; t01 = 0; t11= 0; t21= 0; t12= 0; t22= 0; C12Numerator = 0;
	
		float Sharedace, C12, part1, part2, part3, part4, part5, Gamma1, Gamma2, Gamma3;
	
		/*fARare = number of OTUs with one individual found in A and less than or equal to 10 in B. 
		fBRare = number of OTUs with one individual found in B and less than or equal to 10 in A. 
		arare = number of sequences from A that contain less than 10 sequences. 
		brare = number of sequences from B that contain less than 10 sequences. 
		S12Rare = number of shared OTUs where both of the communities are represented by less than or equal to 10 sequences 
		S12Abund = number of shared OTUs where at least one of the communities is represented by more than 10 sequences 
		S12 = number of shared OTUs in A and B
		This estimator was changed to reflect Caldwell's changes, eliminating the nrare / nrare - 1 */

		for (int i = 0; i < shared[0]->size(); i++) {
			//store in temps to avoid multiple repetitive function calls
			tempA = shared[0]->getAbundance(i);
			tempB = shared[1]->getAbundance(i);
			if ((tempA != 0) && (tempB != 0)) {//they are shared
				S12++;
				//do both A and B have one
				if ((tempA == 1) && (tempB == 1))		{	f11++;	 }
				//is A one and B rare
				if ((tempA == 1) && (tempB <= abund))	{  fARare++; }
				//is B one and A rare
				if ((tempB == 1) && (tempA <= abund))	{  fBRare++; }
			
				if ((tempA <= abund) && (tempB <= abund)) { //shared and both rare
					S12Rare++;
					t10 += tempA;	//Sum Xi
					t01 += tempB;	//Sum Yi
					
					//calculate top of C12
					// YiI(Xi = 1)
					if (tempA == 1) { C12Numerator += tempB; }
					//XiI(Yi = 1)
					if (tempB == 1)	{ C12Numerator += tempA; }
					//-I(Xi=Yi=1)
					if ((tempA == 1) && (tempB == 1)) { C12Numerator--; }
					
					//calculate t11 - Sum of XiYi
					t11 += tempA * tempB;
					//calculate t21  - Sum of Xi(Xi - 1)Yi
					t21 += tempA * (tempA - 1) * tempB;
					//calculate t12  - Sum of Xi(Yi - 1)Yi
					t12 += tempA * (tempB - 1) * tempB;
					//calculate t22  - Sum of Xi(Xi - 1)Yi(Yi - 1)
					t22 += tempA * (tempA - 1) * tempB * (tempB - 1);

				}
				if ((tempA > 10) || (tempB > 10)) {
					S12Abund++;
				}
			}
		}
	
		C12 = 1 - (C12Numerator /(float) t11);
		part1 = S12Rare / (float)C12;
		part2 = 1 / (float)C12;

		//calculate gammas
		Gamma1 = ((S12Rare * t21) / (float)((C12 * t10 * t11)) - 1);
		Gamma2 = ((S12Rare * t12) / (float)((C12 * t01 * t11)) - 1);
		Gamma3 = ((S12Rare / C12) * (S12Rare / C12)) * ( t22 / (float)(t10 * t01 * t11));
		Gamma3 = Gamma3 - ((S12Rare * t11) / (float)(C12 * t01 * t10)) - Gamma1 - Gamma2;	

		if (isnan(Gamma1) || isinf(Gamma1)) { Gamma1 = 0; }
		if (isnan(Gamma2) || isinf(Gamma2)) { Gamma2 = 0; }
		if (isnan(Gamma3) || isinf(Gamma3)) { Gamma3 = 0; }
		if (isnan(part1)  || isinf(part1))  { part1 = 0;  }
		if (isnan(part2)  || isinf(part2))  { part2 = 0;  }
			
		part3 = fARare * Gamma1;
		part4 = fBRare * Gamma2;
		part5 = f11 * Gamma3;

		Sharedace = S12Abund + part1 + (part2 * (part3 + part4 + part5));
		data[0] = Sharedace;
	
		return data;
	}
	catch(exception& e) {
		errorOut(e, "SharedAce", "getValues");
		exit(1);
	}
}
/***********************************************************************/


