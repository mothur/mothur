#ifndef NSEQS_H
#define NSEQS_H

/*
 *  nseqs.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/16/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */



#include "calculator.h"

/***********************************************************************/

class NSeqs : public Calculator {

public:
	NSeqs() : Calculator("nseqs", 1, false) {};
	EstOutput getValues(SAbundVector* rank){
		data.resize(1,0);
		data[0] = (double)rank->getNumSeqs();
		return data;
	}
	
	EstOutput getValues(vector<SharedRAbundVector*> shared) { //return number of sequences in the sharedotus
		
		int numGroups = shared.size();
		data.clear(); data.resize(numGroups,0);

		for (int i = 0; i < shared[0]->getNumBins(); i++) {
			//get bin values and set sharedByAll 
			bool sharedByAll = true;
			for (int j = 0; j < numGroups; j++) {
				if (shared[j]->get(i) == 0) { sharedByAll = false; }
			}
			
			//they are shared
			if (sharedByAll == true) {  for (int j = 0; j < numGroups; j++) {  data[j] += shared[j]->get(i);  } }
		}

		return data;
	}
	string getCitation() { return "http://www.mothur.org/wiki/Nseqs"; }
};

/***********************************************************************/

#endif
