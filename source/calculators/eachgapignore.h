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
 
 
#include "calculator.h"

/**************************************************************************************************/

class eachGapIgnoreTermGapDist : public DistCalc {
	
public:
	eachGapIgnoreTermGapDist() {}
	eachGapIgnoreTermGapDist(const eachGapIgnoreTermGapDist& ddb) {}
	
	double calcDist(Sequence A, Sequence B){
        try {
            int diff = 0;
            int length = 0;
            int start = 0;
            int end = 0;
            bool overlap = false;
            
            string seqA = A.getAligned();
            string seqB = B.getAligned();
            int alignLength = seqA.length();
            
            for(int i=0;i<alignLength;i++){
                if(seqA[i] != '.' && seqB[i] != '.' && seqA[i] != '-' && seqB[i] != '-' ){
                    start = i;
                    overlap = true;
                    break;
                }
            }
            for(int i=alignLength-1;i>=0;i--){
                if(seqA[i] != '.' && seqB[i] != '.' && seqA[i] != '-' && seqB[i] != '-' ){
                    end = i;
                    overlap = true;
                    break;
                }
            }
            
            for(int i=start;i<=end;i++){
                if(seqA[i] == '.' || seqB[i] == '.'){
                    break;
                }
                else if(seqA[i] != '-' || seqB[i] != '-'){
                    if(seqA[i] != seqB[i]){
                        diff++;
                    }
                    length++;
                }
            }
            
            //non-overlapping sequences
            if (!overlap) { length = 0; }
            
            if(length == 0)	{	dist = 1.0000;								}
            else			{	dist = ((double)diff  / (double)length);	}
            
            return dist;
        }
        catch(exception& e) {
            m->errorOut(e, "eachGapIgnoreTermGapDist", "calcDist");
            exit(1);
        }
	}
	
};
/**************************************************************************************************/

#endif


