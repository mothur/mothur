/*
 *  boneh.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/13/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "boneh.h"
#include <math.h>

/***********************************************************************/

//This solves for the value of 'v' using a binary search.
double Boneh::getV(double f1, double n, double rs) {

	if(rs == 0)
		return 0;
	
	double accuracy = .0001;
	double v = 100000.0;
	double step = v/2;
	double ls = v * (1 - pow((1 - f1/(n*v)), n));
	
	while(abs(ls - rs) > accuracy) {
		if(ls > rs)	{	v -= step;	}
		else		{	v += step;	}
		
		ls = v * (1 - pow((1 - f1/(n * v)), n));
		step /= 2;
	}

	return v;
}
	
/***********************************************************************/	
EstOutput Boneh::getValues(SAbundVector* sabund){

	try {
		data.resize(1,0);
		

		bool valid = false;
		double sum = 0;
		double n = (double)sabund->getNumSeqs();
		if(m==0){	m=n;	}
		
		double f1 = (double)sabund->get(1);
		
		for(int i = 1; i < sabund->size(); i++){
			sum += (double)sabund->get(i) * exp(-i);
		}

		if(sabund->get(1) > sum)
			valid = true;
		
		sum = 0;
		if(valid) {
			for(int j = 1; j < sabund->size(); j++){
				sum += sabund->get(j) * pow((1 - (double)j / n), n);
			}
			
			double v = getV(f1, n, sum);
			
			sum = 0;
			for(int j = 1; j < sabund->size(); j++) {
				for (int i = 0; i < sabund->get(j); i++) {
					sum += pow(1 - j / n, n) * (1 - pow(1 - j / n, m));
				}
			}
			sum +=  v * pow(1 - f1/(n*v), n) * (1 - pow(1 - f1/(n*v), m));
		}

		data[0] = sum;
		
		return data;
	}
	catch(exception& e) {
		errorOut(e, "Boneh", "getValues");
		exit(1);
	}
}


/***********************************************************************/
