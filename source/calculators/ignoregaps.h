#ifndef IGNOREGAPS_H
#define IGNOREGAPS_H
/*
 *  ignoregaps.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"

/**************************************************************************************************/

//	this class calculates distances by ignoring all gap characters.  so if seq a has an "A" and seq
//	b has a '-', there is no penalty

class ignoreGaps : public DistCalc {
	
public:
	
    ignoreGaps(double c) : DistCalc(c) {}
	
    vector<double> calcDist(Sequence A, classifierOTU otu, vector<int> cols);
     
    double calcDist(Sequence A, Sequence B); //calc distance between 2 seqeunces
    
private:
    
    vector<int> setStarts(classifierOTU seqA, classifierOTU otu, vector<int> cols);
    vector<int> setEnds(classifierOTU seqA, classifierOTU otu, vector<int> cols);
    
	
};

/**************************************************************************************************/
#endif

