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
 
 
#include "dist.h"

/**************************************************************************************************/

class eachGapIgnoreTermGapDist : public Dist {
	
public:
	
	void calcDist(Sequence A, Sequence B){		
		int diff = 0;
		int length = 0;
		int start = 0;
		int end = 0;
		
		for(int i=0;i<A.getLength();i++){
			if(A.getAligned()[i] == '.' || B.getAligned()[i] == '.' || A.getAligned()[i] == '-' || B.getAligned()[i] == '-'){
			}
			else{
				start = i;
				break;
			}
		}
		for(int i=A.getLength()-1;i>=0;i--){
			if(A.getAligned()[i] == '.' || B.getAligned()[i] == '.' || A.getAligned()[i] == '-' || B.getAligned()[i] == '-'){
			}
			else{
				end = i;
				break;
			}
		}
		
		for(int i=start;i<=end;i++){
			if(A.getAligned()[i] == '-' && B.getAligned()[i] == '-'){}
			else if(A.getAligned()[i] == '.' || B.getAligned()[i] == '.'){
				break;	
			}
			else if(A.getAligned()[i] != '-' && B.getAligned()[i] != '-'){
				if(A.getAligned()[i] != B.getAligned()[i]){
					diff++;
				}
				length++;
			}
			else if(A.getAligned()[i] != '-' || B.getAligned()[i] != '-'){
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


