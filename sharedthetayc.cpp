/*
 *  sharedthetayc.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedthetayc.h"

/***********************************************************************/
EstOutput ThetaYC::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
	try {	
		data.resize(1,0);
		
		int Atotal = 0;
		int Btotal = 0;
		float thetaYC = 0;
		float relA = 0;
		float relB = 0;
		float a = 0;
		float b = 0;
		
		//get the total values we need to calculate the theta denominator sums
		for (int i = 0; i < shared1->size(); i++) {
			//store in temps to avoid multiple repetitive function calls
			Atotal += shared1->getAbundance(i);
			Btotal += shared2->getAbundance(i);
		}
		
		//calculate the theta denominator sums
		for (int j = 0; j < shared1->size(); j++) {
			//store in temps to avoid multiple repetitive function calls
			relA = shared1->getAbundance(j) / (float)Atotal;
			relB = shared2->getAbundance(j) / (float)Btotal;
					
			a += relA * relB;
			b += pow((relA-relB),2);
		}

		thetaYC = a / (float) (b+a);
		
		if (isnan(thetaYC) || isinf(thetaYC)) { thetaYC = 0; }
		
		data[0] = thetaYC;
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ThetaYC class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ThetaYC class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
