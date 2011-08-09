#ifndef EACHGAPDIST_H
#define EACHGAPDIST_H
/*
 *  eachgapdist.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/7/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "dist.h"

/**************************************************************************************************/

class eachGapDist : public Dist {
	
public:
	
	eachGapDist() {}
	eachGapDist(const eachGapDist& ddb) {}
	
	void calcDist(Sequence A, Sequence B){		
		int diff = 0;
		int length = 0;
		int start = 0;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();

		int alignLength = seqA.length();
		
		for(int i=0; i<alignLength; i++){
			if(seqA[i] != '.' || seqB[i] != '.'){
				start = i;
				break;
			}
		}

		for(int i=start;i<alignLength;i++){
			if(seqA[i] == '.' && seqB[i] == '.'){
				break;	
			}
			else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-')){;}
			else{
				if(seqA[i] != seqB[i]){
					diff++;
				}
				length++;
			}
		}
		
		if(length == 0)	{	dist = 1.0000;								}
		else			{	dist = ((double)diff  / (double)length);	}
	}
};

/**************************************************************************************************/

#endif
