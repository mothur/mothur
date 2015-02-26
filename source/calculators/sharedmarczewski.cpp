/*
 *  sharedmarczewski.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedmarczewski.h"

EstOutput SharedMarczewski::getValues(vector<SharedRAbundVector*> vectorShared){
	try {
		SharedRAbundVector* shared1 = vectorShared[0];
		SharedRAbundVector* shared2 = vectorShared[1];
		
		data.resize(1,0);
		
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
		data[0] = (b+c)/(a+b+c);

		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedMarczewski", "getValues");
		exit(1);
	}
}

/***********************************************************************/
