/*
 *  invsimpson.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 8/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "invsimpson.h"
#include "simpson.h"

/***********************************************************************/

EstOutput InvSimpson::getValues(SAbundVector* rank){
	try {
		//vector<double> simpsonData(3,0);
		data.resize(3,0);
		vector<double> simpData(3,0);
		Simpson* simp = new Simpson();
		simpData = simp->getValues(rank);
		
		if(simpData[0] != 0){
			data[0] = 1/simpData[0];
			data[1] = 1/simpData[2];
			data[2] = 1/simpData[1];
		}
		else{
			data.assign(3,1);
		}
		
		delete simp;
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "InvSimpson", "getValues");
		exit(1);
	}
}

/***********************************************************************/
