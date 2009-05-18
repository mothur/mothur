/*
 *  shen.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/18/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "shen.h"
#include <math.h>

	
/***********************************************************************/	
EstOutput Shen::getValues(SAbundVector* rank){

	try {
		data.resize(1,0);
		
		double n = (double)rank->getNumSeqs();
		double f1 = (double)rank->get(1);
		
		double D_rare = 0; //I didn't know what this was. Simply replace the '0' with the appropriate expression.
		double C_rare = 1; //I didn't know what this was. Simply replace the '1' with the appropriate expression.
		double Y_rare = 1; //I didn't know what this was. Simply replace the '1' with the appropriate expression.
		
		double f0 = D_rare/C_rare + f1/C_rare * Y_rare*Y_rare - D_rare;
		
		data[0] = f0 * (1 - pow(1 - f1/n/f0, m));

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
