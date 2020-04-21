//
//  eachgapignore.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/21/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "eachgapignore.h"

/***********************************************************************/

double eachGapIgnoreTermGapDist::calcDist(Sequence A, Sequence B){
    try {
        int diff = 0;
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
        
        //non-overlapping sequences
        if (!overlap) { return 1.0000; }
        
        int maxMinLength = end - start + 1;
        
        for(int i=start;i<alignLength;i++){
            if(seqA[i] == '.' || seqB[i] == '.'){ //reached terminal gaps, so quit
                break;
            }
            else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-')){ maxMinLength--; } //comparing gaps, ignore
            else{
                if(seqA[i] != seqB[i]){
                    diff++;
                }
            }
            
            dist = (double)diff / maxMinLength;
            
            if (dist > cutoff) { return 1.0000; }
        }
        
        if(maxMinLength == 0)    {    dist = 1.0000;                                }
        else            {    dist = ((double)diff  / (double)maxMinLength);        }
        
        return dist;
    }
    catch(exception& e) {
        m->errorOut(e, "eachGapIgnoreTermGapDist", "calcDist");
        exit(1);
    }
}
/***********************************************************************/
vector<double> eachGapIgnoreTermGapDist::calcDist(Sequence A, classifierOTU otu, vector<int> cols){ //this function calcs the distance using only the columns provided
    try {
        vector<double> dists; dists.resize(otu.numSeqs, 0.0);
        
        //if you didn't select columns, use all columns
        if (cols.size() == 0) {
            for (int i = 0; i < otu.otuData.size(); i++) { cols.push_back(i); }
        }
        
        classifierOTU seq(A.getAligned());
        vector<int> starts = setStartsIgnoreTermGap(seq, otu, cols);
        vector<int> ends = setEndsIgnoreTermGap(seq, otu, cols);
        
        int alignLength = cols.size();
        
        for (int h = 0; h < otu.numSeqs; h++) {
            
            if (m->getControl_pressed()) { break; }
            
            if ((starts[h] == -1) && (ends[h] == -1)) { dists[h] = 1.0000; } //no overlap
            else {
                
                if (starts[h] == -1)    { starts[h] = 0;    }
                if (ends[h] == -1)      { ends[h] = 0;      }
                
                int maxMinLength = ends[h] - starts[h] + 1;
                int difference = 0;
                
                for(int i=starts[h];i<alignLength;i++){
                    
                    char seqA = seq.otuData[cols[i]][0];
                    vector<char> otuChars = otu.otuData[cols[i]];
                    char seqB = otuChars[0]; //assume column if identical
                    if (otuChars.size() == otu.numSeqs) {  seqB = otuChars[h]; }
                    
                    if(seqA == '.' || seqB == '.'){ i+=alignLength; } //terminal gaps, break;
                    
                    else if((seqA == '-' && seqB == '-') || (seqA == '-' && seqB == '.') || (seqA == '.' && seqB == '-')){ maxMinLength--; } //comparing gaps, ignore
                    else{
                        if(seqA != seqB){ difference++; }
                    }
                    
                    double distance = 1.0;
                    distance = (double)difference / maxMinLength;
                    
                    if (distance > cutoff) { dists[h] = 1.0000;    i+=alignLength;  } //break;
                }
                
                if(maxMinLength == 0)       {    dists[h] = 1.0000;                               }
                else if (dists[h] == 0.0)   {    dists[h] = (double)difference / maxMinLength;    } //not set
            }
        }
            
        return dists;
     }
     catch(exception& e) {
         m->errorOut(e, "eachGapIgnoreTermGapDist", "calcDist");
         exit(1);
     }
 }

/***********************************************************************/
