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

#include "dist.h"

/**************************************************************************************************/

//	this class calculates distances by ignoring all gap characters.  so if seq a has an "A" and seq
//	b has a '-', there is no penalty

class ignoreGaps : public Dist {
	
public:
	
	ignoreGaps() {}
	
	void calcDist(Sequence A, Sequence B){		
		int diff = 0;
		int length = 0;
		int start = 0;
		bool overlap = false;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();
		int alignLength = seqA.length();
		
		for(int i=0;i<alignLength;i++){
			if(seqA[i] != '.' && seqB[i] != '.'){
				start = i;
				overlap = true;
				break;
			}
		}
		
		for(int i=start; i<alignLength; i++){
			if(seqA[i] == '.' || seqB[i] == '.'){
				break;
			}
			else if((seqA[i] != '-' && seqB[i] != '-')){
				if(seqA[i] != seqB[i]){
					diff++;
				}
				length++;
			}
		}
		
		//non-overlapping sequences
		if (!overlap) { length = 0; }

		if(length == 0)		{	dist = 1.0000;								}
		else				{	dist = ((double)diff  / (double)length);	}
		
	}
	
};

/**************************************************************************************************/
#endif

