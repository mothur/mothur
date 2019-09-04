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


#include "calculator.h"

/**************************************************************************************************/

class eachGapDist : public DistCalc {
	
public:
	
	eachGapDist() {}
    eachGapDist(double c) : DistCalc(c) {}
	
	double calcDist(Sequence A, Sequence B){
		int diff = 0;
        int start = 0;
		
		string seqA = A.getAligned();
		string seqB = B.getAligned();

		int alignLength = (int)seqA.length();
		
		for(int i=0; i<alignLength; i++){
			if(seqA[i] != '.' || seqB[i] != '.'){
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
                if(seqA[i] == '.' && seqB[i] == '.'){ //reached terminal gaps, so quit
                    break;
                }
                else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-')){ maxMinLength--; } //comparing gaps, ignore
                else{
                    if(seqA[i] != seqB[i]){
                        diff++;
                    }
                    //length++;
                }
                
                dist = (double)diff / maxMinLength;
                
                if (dist > cutoff) { return 1.0000; }
            }
            
            if(maxMinLength == 0)	{	dist = 1.0000;								}
            else			{	dist = ((double)diff  / (double)maxMinLength);	    }
            
        }else {
        
            int length = 0;
            
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
            
            if(length == 0)    {    dist = 1.0000;                                }
            else            {    dist = ((double)diff  / (double)length);    }
        }
        
        return dist;
	}
};

/**************************************************************************************************/

#endif
