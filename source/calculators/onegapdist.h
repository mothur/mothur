#ifndef ONEGAPDIST_H
#define ONEGAPDIST_H
/*
 *  onegapdist.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "calculator.h"


/**************************************************************************************************/

class oneGapDist : public DistCalc {
	
public:
	
    oneGapDist(double c) : DistCalc(c) {}
	
    //finds the distance from A to each seq in otu.
    //this function calcs the distance using only the columns provided, if cols is empty, use all
    vector<double> calcDist(Sequence A, classifierOTU otu, vector<int> cols);
    
    double calcDist(Sequence A, Sequence B); //calc distance between 2 seqeunces
    
    string getCitation() { return "http://mothur.org"; }
	
};

/**************************************************************************************************/

#endif
