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
	
	ignoreGaps() {}
    ignoreGaps(double c) : DistCalc(c) {}
	
	double calcDist(Sequence A, Sequence B){
		int diff = 0;
        int start = 0;
		bool overlap = false;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();
		int alignLength = (int)seqA.length();
		
		for(int i=0;i<alignLength;i++){
			if(seqA[i] != '.' && seqB[i] != '.'){
				start = i;
				overlap = true;
				break;
			}
		}
        
        if (setCutoff) {
            
            int end = 0;
            for(int i=alignLength-1;i>=0;i--){
                if(seqA[i] != '.' && seqB[i] != '.'){
                    end = i;
                    overlap = true;
                    break;
                }
            }
            
            //non-overlapping sequences
            if (!overlap) { return 1.0000; }
            
            int maxMinLength = end - start + 1;
            
            for(int i=start; i<alignLength; i++){
                if(seqA[i] == '.' || seqB[i] == '.'){
                    break;
                }
                else if((seqA[i] != '-' && seqB[i] != '-')){
                    if(seqA[i] != seqB[i]){
                        diff++;
                    }
                }else {  maxMinLength--; }
                
                dist = (double)diff / maxMinLength;
                
                if (dist > cutoff) { return 1.0000; }
            }
            
            if(maxMinLength == 0)	{	dist = 1.0000;								    }
            else				    {	dist = ((double)diff  / (double)maxMinLength);	}
            
        }else { //for lt options where we need the dists above the cutoff
            int length = 0;
            
            //non-overlapping sequences
            if (!overlap) { return 1.0000; }
            
            for(int i=start; i<alignLength; i++){
                if(seqA[i] == '.' || seqB[i] == '.'){ break; }
                
                else if((seqA[i] != '-' && seqB[i] != '-')){
                    
                    if(seqA[i] != seqB[i]){ diff++; }
                    length++;
                }
            }
            
            if(length == 0)        {    dist = 1.0000;                                }
            else                {    dist = ((double)diff  / (double)length);    }
        }
        
        return dist;
		
	}
	
};

/**************************************************************************************************/
#endif

