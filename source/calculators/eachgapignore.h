#ifndef EACHGAPIGNORE_H
#define EACHGAPIGNORE_H
/*
 *  eachgapignore.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
 
#include "calculator.h"

/**************************************************************************************************/

class eachGapIgnoreTermGapDist : public DistCalc {
	
public:
    
    eachGapIgnoreTermGapDist(double c) : DistCalc(c) {}
	
    double calcDist(Sequence A, Sequence B);
	
    vector<double> calcDist(Sequence A, classifierOTU otu, vector<int> cols);

	
};
/**************************************************************************************************/

#endif


