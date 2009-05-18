/*
 *  boneh.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "boneh.h"
#include <math.h>

/***********************************************************************/
//This solves for the value of 'v' using a binary search.
double Boneh::getV(double f1, double n, double rs) {

	//cout << "f1 = " << f1 << "\nn = " << n << "\nrs = " << rs << "\n\n";
	
	if(rs == 0)
		return 0;
	
	double accuracy = .0001;
	double v = 100000.0;
	double step = v/2;
	double ls = v * (1 - pow((1 - f1/(n*v)), n));
	
	//cout << "ls = " << ls << "\n";
	
	while(abs(ls - rs) > accuracy) {
		if(ls > rs)
			v -= step;
		else
			v += step;
		
		ls = v* (1 - pow((1 - f1/(n*v)), n));
		step /= 2;
		
		//cout << "ls = " << ls << "\n";
	}
	
	return v;
}
	
/***********************************************************************/	
EstOutput Boneh::getValues(SAbundVector* rank){

	try {
		data.resize(1,0);
		
		bool valid = false;
		double sum = 0;
		double n = (double)rank->getNumSeqs();
		double f1 = (double)rank->get(1);
		
		for(int i = 1; i < rank->size(); i++)
			sum += (double)rank->get(i) * exp(-i);
		
		if(rank->get(1) > sum)
			valid = true;
		
		sum = 0;
		if(valid) {
			for(int j = 1; j < rank->size(); j++)
				sum += rank->get(j) * pow((1 - (double)j / n), n);
			
			double v = getV(f1, n, sum);
			
			//cout << "v = " << v << "\n";
			
			
			sum = 0;
			for(int j = 1; j < rank->size(); j++) {
				double Xi = 0; //I didn't know what this was, simply replace the 0
							   //with the appropriate expression for the boneh calculator
							   //to work.
				sum += pow(1 - Xi / n, n) * (1 - pow(1 - Xi / n, m)) + v * pow(1 - f1/(n*v), n) * (1 - pow(1 - f1/(n*v), m));
			}
		}

		data[0] = sum;
		
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Coverage class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Coverage class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
};


/***********************************************************************/
