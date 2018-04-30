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
	
	oneGapDist() {}
	
	double calcDist(Sequence A, Sequence B){
		
		int difference = 0;
		int minLength = 0;
		int openGapA = 0;
		int openGapB = 0;
		int start = 0;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();
		int alignLength = seqA.length();
		
		for(int i=0;i<alignLength;i++){
			if((seqA[i] != '.' || seqB[i] != '.')){
				start = i;
				break;
			}
		}

		for(int i=start;i<alignLength;i++){
			if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '.' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.')){	;	}
			else if(seqA[i] == '.' && seqB[i] == '.'){
				break;
			}
			else if(seqB[i] != '-' && (seqA[i] == '-' || seqA[i] == '.')){
				if(openGapA == 0){
					difference++;
					minLength++;
					openGapA = 1;
					openGapB = 0;
				}
			}
			else if(seqA[i] != '-' && (seqB[i] == '-' || seqB[i] == '.')){
				if(openGapB == 0){
					difference++;
					minLength++;
					openGapA = 0;
					openGapB = 1;
				}
			}
			else if(seqA[i] != '-' && seqB[i] != '-'){
				if(seqA[i] != seqB[i]){
					difference++;
					minLength++;
					openGapA = 0;
					openGapB = 0;
				}
				else{
					minLength++;
					openGapA = 0;
					openGapB = 0;
				}
			}
		}
	
		if(minLength == 0)	{	dist = 1.0000;							}
		else				{	dist = (double)difference / minLength;	}
        return dist;
	}
	
};

/**************************************************************************************************/

#endif
