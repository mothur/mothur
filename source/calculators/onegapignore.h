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
	
	double calcDist(Sequence A, Sequence B){
		
		int difference = 0;
		bool openGapA = false;
		bool openGapB = false;
		int start = 0;
		int end = 0;
		bool overlap = false;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();
		int alignLength = (int)seqA.length();

		// this assumes that sequences start and end with '.'s instead of'-'s.
		for(int i=0;i<alignLength;i++){
			if(seqA[i] != '.' && seqB[i] != '.' && seqA[i] != '-' && seqB[i] != '-' ){ //skip leading gaps
				start = i;

				overlap = true;
				break;
			}
		}
		for(int i=alignLength-1;i>=0;i--){
			if(seqA[i] != '.' && seqB[i] != '.' && seqA[i] != '-' && seqB[i] != '-' ){ //ignore terminal gaps
				end = i;

				overlap = true;
				break;
			}
		}
		
        //non-overlapping sequences
        if (!overlap) { return 1.0000; }
        
        if (setCutoff) {
            int maxMinLength = end - start;
            
            for(int i=start;i<=end;i++){
                if(seqA[i] == '-' && seqB[i] == '-'){	maxMinLength--;	} //comparing gaps, ignore
                
                else if(seqB[i] != '-' && seqA[i] == '-'){ //seqB is a base, seqA is a gap
                    if(!openGapA){
                        difference++;
                        openGapA = true;
                        openGapB = false;
                    }else { maxMinLength--; }
                }
                else if(seqA[i] != '-' && seqB[i] == '-'){ //seqA is a base, seqB is a gap
                    if(!openGapB){
                        difference++;
                        openGapA = false;
                        openGapB = true;
                    }else { maxMinLength--; }
                }
                else if(seqA[i] != '-' && seqB[i] != '-'){ //both bases
                    openGapA = false;
                    openGapB = false;
                    
                    //no match
                    if(seqA[i] != seqB[i]){ difference++; }
                }
                
                dist = (double)difference / maxMinLength;
                
                if (dist > cutoff) { return 1.0000; }
            }
            
            if(maxMinLength == 0)	{	dist = 1.0000;							    }
            else				    {	dist = (double)difference / maxMinLength;	}
            
        }else {
            int minLength = 0;
            
            for(int i=start;i<=end;i++){
                if(seqA[i] == '-' && seqB[i] == '-'){    ;    }
                else if(seqB[i] != '-' && seqA[i] == '-'){
                    if(!openGapA){
                        difference++;
                        minLength++;
                        openGapA = true;
                        openGapB = false;
                    }
                }
                else if(seqA[i] != '-' && seqB[i] == '-'){
                    if(!openGapB){
                        difference++;
                        minLength++;
                        openGapA = false;
                        openGapB = true;
                    }
                }
                else if(seqA[i] != '-' && seqB[i] != '-'){
                    if(seqA[i] != seqB[i]){
                        difference++;
                    }
                    minLength++;
                    openGapA = false;
                    openGapB = false;
                }
            }
            
            if(minLength == 0)  {    dist = 1.0000;                            }
            else                {    dist = (double)difference / minLength;    }
        }
        return dist;
	}

};

/**************************************************************************************************/

#endif

