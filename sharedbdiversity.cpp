/*
 *  bdiversity.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 3/13/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "sharedbdiversity.h"

/***********************************************************************/

double SharedBDiversity::whitt(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	double nz1 = (double)shared1->numNZ();
	double nz2 = (double)shared2->numNZ();
	double sum = nz1;
	for(int i = 1; i < shared1->size(); i++)
	{
		if(shared1->get(i).abundance == 0 && shared2->get(i).abundance != 0)
			sum++;
			}
	return 2*sum/(nz1+nz2)-1;
}

/***********************************************************************/

double SharedBDiversity::ms(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	double a = 0;
	double b = 0;
	double c = 0;
	for(int i = 1; i < shared1->size(); i++)
	{
		int abund1 = shared1->get(i).abundance;
		int abund2 = shared2->get(i).abundance;
		
		if(abund1 > 0 && abund2 > 0)
			a++;
		else if(abund1 > 0 && abund2 == 0)
			b++;
		else if(abund1 == 0 && abund2 > 0)
			c++;
	}
	return (b+c)/(a+b+c);
}

/***********************************************************************/

double SharedBDiversity::sor(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	double sum = 0;
	double asum = 0;
	double bsum = 0;
	for(int i = 1; i < shared1->size(); i++)
	{
		int abund1 = shared1->get(i).abundance;
		int abund2 = shared2->get(i).abundance;
		
		asum += abund1;
		bsum += abund2;
		if(abund1 >= abund2)
			sum += abund2;
		else
			sum += abund1;
	}
	return 2*sum/(asum+bsum);
}

/***********************************************************************/

double SharedBDiversity::mor(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	double multSum = 0;
	double powSum1 = 0;
	double powSum2 = 0;
	double ind1 = 0;
	double ind2 = 0;
	for(int i = 1; i < shared1->size(); i++)
	{
		double abund1 = (double)shared1->get(i).abundance;
		double abund2 = (double)shared2->get(i).abundance;
		multSum += abund1*abund2;
		powSum1 += pow(abund1, 2);
		powSum2 += pow(abund2, 2);
		ind1 += abund1;
		ind2 += abund2;
	}
	return 2*multSum / ((powSum1/pow(ind1, 2) + powSum2/pow(ind2, 2)) * (ind1*ind2));
}
/***********************************************************************/
	
EstOutput SharedBDiversity::getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2){
	try {
		data.resize(4,0);
		
		double whitt1 = whitt(shared1, shared2);
		double jac1 = 1-ms(shared1, shared2);
		double sor1 = sor(shared1, shared2);
		double mor1 = mor(shared1, shared2);
		/*cout << "Whittaker's Measure: " << whitt1 << "\n";
		cout << "Marczewski-Steinhaus distance: " << jac1 << "\n";
		cout << "Sorensen index: " << sor1 << "\n";
		cout << "Morisita-Horn index: " << mor1 << "\n";*/
		
		data[0] = whitt1;
		data[1] = jac1;
		data[2] = sor1;
		data[3] = mor1;
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		if (isnan(data[1]) || isinf(data[1])) { data[1] = 0; }
		if (isnan(data[2]) || isinf(data[0])) { data[2] = 0; }
		if (isnan(data[3]) || isinf(data[1])) { data[3] = 0; }
		
		return data;
	}
		
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the NPShannon class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the NPShannon class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/