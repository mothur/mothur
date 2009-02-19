#ifndef UVEST_H
#define UVEST_H
/*
 *  uvest.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
/* This class implements the UVEst estimator on two groups. 
It is used by sharedJAbund and SharedSorensonAbund. */
 
 
using namespace std;

#include "mothur.h"
#include "sharedrabundvector.h"

typedef vector<double> EstOutput;

/***********************************************************************/
class UVEst {
	public:
		EstOutput getUVest(SharedRAbundVector* shared1, SharedRAbundVector* shared2);		
};
/***********************************************************************/

#endif
