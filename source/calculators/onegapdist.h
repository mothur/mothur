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
	
    oneGapDist() {} //sets cutoff to 1.0
    oneGapDist(double c) : DistCalc(c) {}
	
	double calcDist(Sequence A, Sequence B){
		
		int difference = 0;
		bool openGapA = false;
		bool openGapB = false;
        int start = 0;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();
		int alignLength = seqA.length();
		
		for(int i=0;i<alignLength;i++){
			if((seqA[i] != '.' || seqB[i] != '.')){ //one of you is not a terminal gap
				start = i;
				break;
			}
		}
        
        if (setCutoff) {
            
            int end = 0;
            for(int i=alignLength-1;i>=0;i--){
                if((seqA[i] != '.' || seqB[i] != '.')){ //one of you is not a terminal gap
                    end = i;
                    break;
                }
            }
            
            int maxMinLength = end - start + 1;
            
            for(int i=start;i<alignLength;i++){
                
                //comparing gaps, ignore
                if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '.' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.')){	maxMinLength--;	}
                //trailing gaps, quit we already calculated all the diffs
                else if(seqA[i] == '.' && seqB[i] == '.'){ break; }
                
                else if(seqB[i] != '-' && (seqA[i] == '-' || seqA[i] == '.')){ //seqB is a base, seqA is a gap
                    if(!openGapA){
                        difference++;
                        openGapA = true;
                        openGapB = false;
                    }else { maxMinLength--; }
                }
                else if(seqA[i] != '-' && (seqB[i] == '-' || seqB[i] == '.')){ //seqA is a base, seqB is a gap
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
            
            
            if(maxMinLength == 0)    {    dist = 1.0000;                          }
            else                {    dist = (double)difference / maxMinLength;    }
        }else {
            int minLength = 0;
            
            for(int i=start;i<alignLength;i++){
                if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '.' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.')){    ;    }
                else if(seqA[i] == '.' && seqB[i] == '.'){
                    break;
                }
                else if(seqB[i] != '-' && (seqA[i] == '-' || seqA[i] == '.')){
                    if(!openGapA){
                        difference++;
                        openGapA = true;
                        openGapB = false;
                        minLength++;
                    }
                }
                else if(seqA[i] != '-' && (seqB[i] == '-' || seqB[i] == '.')){
                    if(!openGapB){
                        difference++;
                        openGapA = false;
                        openGapB = true;
                        minLength++;
                    }
                }
                else if(seqA[i] != '-' && seqB[i] != '-'){
                    minLength++;
                    openGapA = false;
                    openGapB = false;
                    if(seqA[i] != seqB[i]){
                        difference++; }
                }
            }
            
            if(minLength == 0)    {    dist = 1.0000;                            }
            else                {    dist = (double)difference / minLength;    }
        }
        return dist;
	}
	
};

/**************************************************************************************************/

#endif
