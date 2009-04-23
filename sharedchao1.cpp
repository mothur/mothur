/*
 *  sharedchao1.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedchao1.h"

/***********************************************************************/
EstOutput SharedChao1::getValues(SharedRAbundVector* sharedA, SharedRAbundVector* sharedB){
	try {
		data.resize(1,0);
		
		int f11, f1A, f2A, f1B, f2B, S12, tempA, tempB;
		f11 = 0; f1A = 0; f2A = 0; f1B = 0; f2B = 0; S12 = 0;
		float ChaoAB, part1, part2, part3;

	/*	f11 = number of shared OTUs with one observed individual in A and B 
		f1A, f2A = number of shared OTUs with one or two individuals observed in A 
		f1B, f2B = number of shared OTUs with one or two individuals observed in B 
		S12 = number of shared OTUs in A and B  */

		//loop through vectors calculating the f11, f1A, f2A, f1B, f2B, S12 values
		for (int i = 0; i< sharedA->size(); i++) {
			tempA = sharedA->getAbundance(i); //store in temps to avoid calling getAbundance multiple times
			tempB = sharedB->getAbundance(i);
			if ((tempA != 0) && (tempB != 0)) {//they are shared
				S12++; //they are shared
				//do both A and B have one
				if ((tempA == 1) && (tempB == 1)) { f11++; }
				
				//does A have one or two
				if (tempA == 1)			{ f1A++; }
				else if (tempA == 2)	{ f2A++; }
				
				//does B have one or two
				if (tempB == 1)			{ f1B++; }
				else if (tempB == 2)	{ f2B++; }

			}
		}
		
		//checks for divide by zero error
		if ((f2A == 0) || (f2B == 0)) {
			part1 = ((float)(f1A*f1B)/(float)(4*(f2A+1)*(f2B+1)));
			part2 = ((float)(f1A*(f1A-1))/(float)(2*f2A+2));
			part3 = ((float)(f1B*(f1B-1))/(float)(2*f2B+2));
		}else {
			part1 = ((float)(f1A*f1B)/(float)(4*f2A*f2B));
			part2 = ((float)(f1A*f1A)/(float)(2*f2A));
			part3 = ((float)(f1B*f1B)/(float)(2*f2B));
		}
		
		ChaoAB = (float)S12 + (float)(f11*part1) + part2 + part3;
		data[0] = ChaoAB;
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedChao1 class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedChao1 class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	

}

/***********************************************************************/
