//
//  onegapdist.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "onegapdist.h"

/***********************************************************************/

double oneGapDist::calcDist(Sequence A, Sequence B){
    try {
        int difference = 0;
        bool openGapA = false;
        bool openGapB = false;
        
        string seqA = A.getAligned();
        string seqB = B.getAligned();
        int alignLength = seqA.length();
        
        int start = setStart(seqA, seqB);
        int end = setEnd(seqA, seqB);
        
        int maxMinLength = end - start + 1;
        
        for(int i=start;i<alignLength;i++){
            
            //comparing gaps, ignore
            if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '.' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.')){    maxMinLength--;    }
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
        
        return dist;
    }
    catch(exception& e) {
        m->errorOut(e,  "oneGapDist", "calcDist");
        exit(1);
    }
}
/***********************************************************************/

vector<double> oneGapDist::calcDist(Sequence A, classifierOTU otu, vector<int> cols){ //this function calcs the distance using only the columns provided
    try {
        vector<double> dists; dists.resize(otu.numSeqs, 0.0);
        
        //if you didn't select columns, use all columns
        if (cols.size() == 0) {
            for (int i = 0; i < otu.otuData.size(); i++) { cols.push_back(i); }
        }
        
        classifierOTU seq(A.getAligned());
        vector<int> starts = setStarts(seq, otu, cols);
        vector<int> ends = setEnds(seq, otu, cols);
        
        int alignLength = cols.size();

        for (int h = 0; h < otu.numSeqs; h++) {
            
            if (m->getControl_pressed()) { break; }
            
            int maxMinLength = ends[h] - starts[h] + 1;
            int difference = 0;
            bool openGapA = false;
            bool openGapB = false;
            
            for(int i=starts[h];i<alignLength;i++){
                
                char seqA = seq.otuData[cols[i]][0];
                vector<char> otuChars = otu.otuData[cols[i]];
                char seqB = otuChars[0]; //assume column if identical
                if (otuChars.size() == otu.numSeqs) {  seqB = otuChars[h]; }
                
                
                if((seqA == '-' && seqB == '-') || (seqA == '.' && seqB == '-') || (seqA == '-' && seqB == '.')){    maxMinLength--;    }
                
                //trailing gaps, quit we already calculated all the diffs
                else if(seqA == '.' && seqB == '.'){ i+=alignLength; } //break;
                
                else if(seqB != '-' && (seqA == '-' || seqA == '.')){ //seqB is a base, seqA is a gap
                    if(!openGapA){
                        difference++;
                        openGapA = true;
                        openGapB = false;
                    }else { maxMinLength--; }
                }
                
                else if(seqA != '-' && (seqB == '-' || seqB == '.')){ //seqA is a base, seqB is a gap
                    if(!openGapB){
                        difference++;
                        openGapA = false;
                        openGapB = true;
                    }else { maxMinLength--; }
                }
                else if(seqA != '-' && seqB != '-'){ //both bases
                    openGapA = false;
                    openGapB = false;
                    
                    //no match
                    if(seqA != seqB){ difference++; }
                }
                
                double distance = 1.0;
                distance = (double)difference / maxMinLength;
                
                if (distance > cutoff) { dists[h] = 1.0000;    i+=alignLength;  } //break;
            }
            
            if(maxMinLength == 0)       {    dists[h] = 1.0000;                               }
            else if (dists[h] == 0.0)   {    dists[h] = (double)difference / maxMinLength;    } //not set
        }
    
        return dists;
    }
    catch(exception& e) {
        m->errorOut(e, "oneGapDist", "calcDist");
        exit(1);
    }
}
/***********************************************************************/



