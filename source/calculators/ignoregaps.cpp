//
//  ignoregaps.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/21/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "ignoregaps.h"

/***********************************************************************/
double ignoreGaps::calcDist(Sequence A, Sequence B){
    try {
        int diff = 0;
        int start = 0;
        int end = 0;
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
        
        if(maxMinLength == 0)    {    dist = 1.0000;                                    }
        else                    {    dist = ((double)diff  / (double)maxMinLength);    }
        
        return dist;
    }
    catch(exception& e) {
        m->errorOut(e,  "ignoreGaps", "calcDist");
        exit(1);
    }
}
/***********************************************************************/

vector<double> ignoreGaps::calcDist(Sequence A, classifierOTU otu, vector<int> cols){ //this function calcs the distance using only the columns provided
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
            
            if ((starts[h] == -1) && (ends[h] == -1)) { dists[h] = 1.0000; } //no overlap
            else {
                
                if (starts[h] == -1)    { starts[h] = 0;    }
                if (ends[h] == -1)      { ends[h] = 0;      }
                
                int maxMinLength = ends[h] - starts[h] + 1;
                
                int difference = 0;
                bool openGapA = false;
                bool openGapB = false;
                
                for(int i=starts[h];i<alignLength;i++){
                    
                    char seqA = seq.otuData[cols[i]][0];
                    vector<char> otuChars = otu.otuData[cols[i]];
                    char seqB = otuChars[0]; //assume column if identical
                    if (otuChars.size() == otu.numSeqs) {  seqB = otuChars[h]; }
                    
                    if(seqA == '.' || seqB == '.'){ i+=alignLength; }
                    
                    else if((seqA != '-' && seqB != '-')){
                        
                        if(seqA != seqB){ difference++; }
                        
                    }else {  maxMinLength--; }
                    
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
        m->errorOut(e, "ignoreGaps", "calcDist");
        exit(1);
    }
}
/***********************************************************************/
vector<int> ignoreGaps::setStarts(classifierOTU seqA, classifierOTU otu, vector<int> cols){
    try {
        vector<int> starts; starts.resize(otu.numSeqs, -1);
        
        int alignLength = cols.size();
        
        int seqAStart = 0;
        for(int i=0;i<alignLength;i++){ //for each column we want to include
            if (seqA.otuData[cols[i]][0] != '.') {
                seqAStart = i; break;
            }
        }
        
        //set start positions
        int numset = 0;
        for(int i=seqAStart;i<alignLength;i++){ //start can't be before seqAStart because of the &&
            
            if(numset == otu.numSeqs) { break; }
            
            vector<char> thisColumn = otu.otuData[cols[i]];
            if (thisColumn.size() != otu.numSeqs) { //all seqs at this spot are identical
                
                char thisChar = thisColumn[0];
                
                if (thisChar == '.') { } //every seq in otu is a '.' or '-' at this location, move to next column
                else { //this is a base in all locations, you are done
                    for (int k = 0; k < starts.size(); k++) {
                        if ((starts[k] == -1) && (seqA.otuData[cols[i]][0] != '.')) { starts[k] = i; numset++; } //any unset starts are set to this location
                    }
                    break;
                }
            }else{
                for(int j=0;j<otu.numSeqs;j++){ //for each reference
                    if((thisColumn[j] != '.') && (starts[j] == -1) && (seqA.otuData[cols[i]][0] != '.')){ //seq j hasn't set the start value and its a base
                        starts[j] = i; numset++;
                    }
                }
            }
        }
        
        return starts;
    }
    catch(exception& e) {
        m->errorOut(e, "ignoreGaps", "setStarts");
        exit(1);
    }
}
/***********************************************************************/

vector<int> ignoreGaps::setEnds(classifierOTU seqA, classifierOTU otu, vector<int> cols){
    try {
        vector<int> ends; ends.resize(otu.numSeqs, -1);
        
        int alignLength = cols.size();
        
        int seqAEnd = 0;
        for(int i=alignLength-1;i>=0;i--){//for each column we want to include
            if (seqA.otuData[cols[i]][0] != '.') {
                seqAEnd = i; break;
            }
        }
        
        //set start positions
        int numset = 0;
        for(int i=seqAEnd;i>=0;i--){ //for each column we want to include
            
            if(numset == otu.numSeqs) { break; }
            
            vector<char> thisColumn = otu.otuData[cols[i]];
            if (thisColumn.size() != otu.numSeqs) { //all seqs at this spot are identical
                
                char thisChar = thisColumn[0];
                
                if (thisChar == '.'){ } //every seq in otu is a '.' at this location, move to next column
                else { //this is a base in all locations, you are done
                    for (int k = 0; k < ends.size(); k++) {
                        if ((ends[k] == -1) && (seqA.otuData[cols[i]][0] != '.')){ ends[k] = i; numset++; } //any unset starts are set to this location
                    }
                    break;
                }
            }else{
                for(int j=0;j<otu.numSeqs;j++){ //for each reference
                    if((thisColumn[j] != '.') && (ends[j] == -1) && (seqA.otuData[cols[i]][0] != '.')) { //seq j hasn't set the start value and its a base
                        ends[j] = i; numset++;
                    }
                }
            }
        }
        
        return ends;
    }
    catch(exception& e) {
        m->errorOut(e, "ignoreGaps", "setEnds");
        exit(1);
    }
}
/***********************************************************************/
