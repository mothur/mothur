#ifndef EACHGAPDISTIGNORENS_H
#define EACHGAPDISTIGNORENS_H

/*
 *  eachgapdistignorens.h
 *  Mothur
 *
 *  Created by Pat Schloss on 4/20/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "calculator.h"

/**************************************************************************************************/

class eachGapDistIgnoreNs : public DistCalc {
	
public:
	double calcDist(Sequence A, Sequence B){
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
			else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-') || seqA[i] == 'N' || seqB[i] == 'N'){;}
			else{
				if(seqA[i] != seqB[i]){
					diff++;
				}
				length++;
			}
		}
		
		if(length == 0)	{	dist = 1.0000;								}
		else			{	dist = ((double)diff  / (double)length);	}
        
        return dist;
	}
};

/**************************************************************************************************/

#endif
