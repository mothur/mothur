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
    
    eachGapDistIgnoreNs() {}
    eachGapDistIgnoreNs(double c) : DistCalc(c) {}
	
public:
	double calcDist(Sequence A, Sequence B){
		int diff = 0;
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
		
        if (setCutoff) {
            
            int end = 0;
            
            for(int i=alignLength-1;i>=0;i--){
                if(seqA[i] != '.' || seqB[i] != '.'){ //ignore terminal gaps
                    end = i;
                    break;
                }
            }
            
            int maxMinLength = end - start + 1;
            
            for(int i=start;i<alignLength;i++){
                if(seqA[i] == '.' && seqB[i] == '.'){ //reached terminal gaps, so quit
                    break;
                }
                else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-') || seqA[i] == 'N' || seqB[i] == 'N')){ maxMinLength--; } //comparing gaps, ignore
                else{
                    if(seqA[i] != seqB[i]){ diff++; }
                }
                
                dist = (double)diff / maxMinLength;
                
                if (dist > cutoff) { return 1.0000; }
            }
            
            if(maxMinLength == 0)    {    dist = 1.0000;                                }
            else            {    dist = ((double)diff  / (double)maxMinLength);        }
            
        }else {
            int length = 0;
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
        }
        
        return dist;
	}
};

/**************************************************************************************************/

#endif
