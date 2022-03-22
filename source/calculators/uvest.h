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
 

#include "mothurout.h"
#include "calculator.h"

typedef vector<double> EstOutput;

/***********************************************************************/
class UVEst {
	public:
		UVEst() { m = MothurOut::getInstance(); }
        ~UVEst() = default;
    
		EstOutput getUVest(vector<SharedRAbundVector*>);
    
	private:
		MothurOut* m;
};
/***********************************************************************/

#endif
