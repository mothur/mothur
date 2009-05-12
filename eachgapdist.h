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
	void calcDist(Sequence A, Sequence B){		
		int diff = 0;
		int length = 0;
		int start = 0;
		
		for(int i=0; i<A.getLength(); i++){
			if(A.getUnaligned()[i] != '.' || B.getUnaligned()[i] != '.'){
				start = i;
				break;
			}
		}
		
		for(int i=start;i<A.getLength();i++){
			if((A.getUnaligned()[i] == '-' || B.getUnaligned()[i] == '-') && (A.getUnaligned()[i] == '.' || B.getUnaligned()[i] == '.')){}
			else if(A.getUnaligned()[i] == '.' && B.getUnaligned()[i] == '.'){
				break;	
			}
			else if(A.getUnaligned()[i] != '-' && B.getUnaligned()[i] != '-'){
				if(A.getUnaligned()[i] != B.getUnaligned()[i]){
					diff++;
				}
				length++;
			}
			else if(A.getUnaligned()[i] != '-' || B.getUnaligned()[i] != '-'){
				diff++;
				length++;
			}
		}
		
		if(length == 0)	{	dist = 1.0000;								}
		else			{	dist = ((double)diff  / (double)length);	}
		
	}
};

/**************************************************************************************************/

#endif