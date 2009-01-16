/*
 *  rarecalc.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarecalc.h"

/***********************************************************************/

EstOutput RareCalc::getValues(int n){
	try {
		EstOutput output(3,0);

		double richness = (double)numBins;
		double varS = 0.0000;
	
		double varTerm1 = 0.0000;
		double varTerm2 = 0.0000;
	
		double rSummation = 0;
		for(int i=0;i<numBins;i++){
			int N_ni = numSeqs - bins->get(i);
			rSummation += (bMatrix[N_ni][n]);
		
			varTerm1 += (bMatrix[N_ni][n] * (1.0 - bMatrix[N_ni][n] / bMatrix[numSeqs][n]));

			for(int j=i+1;j<numBins;j++){
				varTerm2 += ( bMatrix[N_ni-bins->get(j)][n] - bMatrix[N_ni][n] * bMatrix[numSeqs-bins->get(j)][n] / bMatrix[numSeqs][n]);
			}
		
		}
		richness -= (rSummation / bMatrix[numSeqs][n]);
		varS = (varTerm1 + 2 * varTerm2) / bMatrix[numSeqs][n];
		float sd = pow(varS, 0.5);
	

		output[0] = richness;
		output[1] = richness - 1.96 * sd;	
		output[2] = richness + 1.96 * sd;
	
		return output;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareCalc class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareCalc class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
