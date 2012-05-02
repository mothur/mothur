/*
 *  sharedjabund.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedjabund.h"

/***********************************************************************/

EstOutput JAbund::getValues(vector<SharedRAbundVector*> shared) {
	try {
		EstOutput UVest;
		UVest.resize(2,0);
		data.resize(1,0);
		
		UVest = uv->getUVest(shared);
		
		//UVest[0] is Uest UVest[1] is Vest
		data[0] = 1.0-(UVest[0] * UVest[1]) / ((float)(UVest[0] + UVest[1] - (UVest[0] * UVest[1])));
		if(data[0] > 1){data[0] = 0;	}
		if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "JAbund", "getValues");
		exit(1);
	}
}

/***********************************************************************/
