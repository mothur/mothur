#ifndef ONEIGNOREGAPS_H
#define ONEIGNOREGAPS_H
/*
 *  onegapignore.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "calculator.h"

/**************************************************************************************************/

class oneGapIgnoreTermGapDist : public DistCalc {
	
public:
	
	oneGapIgnoreTermGapDist() {}
    oneGapIgnoreTermGapDist(double c) : DistCalc(c) {}
	
    //finds the distance from A to each seq in otu.
      //this function calcs the distance using only the columns provided, if cols is empty, use all
    vector<double> calcDist(Sequence A, classifierOTU otu, vector<int> cols);
      
    double calcDist(Sequence A, Sequence B); //calc distance between 2 seqeunces
    
private:
       
    
       vector<int> setStarts(classifierOTU seqA, classifierOTU otu, vector<int> cols);
       vector<int> setEnds(classifierOTU seqA, classifierOTU otu, vector<int> cols);
       

};

/**************************************************************************************************/

#endif

