//
//  onegapignore.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/20/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "onegapignore.h"

/***********************************************************************/

double oneGapIgnoreTermGapDist::calcDist(Sequence A, Sequence B){
    
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
            if(seqA[i] == '-' && seqB[i] == '-'){    maxMinLength--;    } //comparing gaps, ignore
            
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
        
        if(maxMinLength == 0)    {    dist = 1.0000;                                }
        else                    {    dist = (double)difference / maxMinLength;    }
        
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
/***********************************************************************/
vector<double> oneGapIgnoreTermGapDist::calcDist(Sequence A, classifierOTU otu, vector<int> cols){ //this function calcs the distance using only the columns provided
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
        
        if (setCutoff) {
            
            for (int h = 0; h < otu.numSeqs; h++) {
                
                if (m->getControl_pressed()) { break; }
                
                if ((starts[h] == -1) && (ends[h] == -1)) { dists[h] = 1.0000; } //no overlap
                else {
                    
                    if (starts[h] == -1)    { starts[h] = 0;    }
                    if (ends[h] == -1)      { ends[h] = 0;      }
                    
                    int maxMinLength = ends[h] - starts[h];
                    int difference = 0;
                    bool openGapA = false;
                    bool openGapB = false;
                    
                    for(int i=starts[h];i<=ends[h];i++){
                        
                        char seqA = seq.otuData[cols[i]][0];
                        vector<char> otuChars = otu.otuData[cols[i]];
                        char seqB = otuChars[0]; //assume column if identical
                        if (otuChars.size() == otu.numSeqs) {  seqB = otuChars[h]; }
                        
                        
                        if(seqA == '-' && seqB == '-'){    maxMinLength--;    } //comparing gaps, ignore
                        
                        else if(seqB != '-' && seqA == '-'){ //seqB is a base, seqA is a gap
                            if(!openGapA){
                                difference++;
                                openGapA = true;
                                openGapB = false;
                            }else { maxMinLength--; }
                        }
                        else if(seqA != '-' && seqB == '-'){ //seqA is a base, seqB is a gap
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
            }
            
        }else{
            for (int h = 0; h < otu.numSeqs; h++) {
                
                if (m->getControl_pressed()) { break; }
                
                int minLength = 0;
                int difference = 0;
                bool openGapA = false;
                bool openGapB = false;
                
                for(int i=starts[h];i<ends[h];i++){
                    
                    char seqA = seq.otuData[cols[i]][0];
                    vector<char> otuChars = otu.otuData[cols[i]];
                    char seqB = otuChars[0]; //assume column if identical
                    if (otuChars.size() == otu.numSeqs) {  seqB = otuChars[h]; }
                    
                    if(seqA == '-' && seqB == '-'){    ;    }
                    else if(seqB != '-' && seqA == '-'){
                        if(!openGapA){
                            difference++;
                            minLength++;
                            openGapA = true;
                            openGapB = false;
                        }
                    }
                    else if(seqA != '-' && seqB == '-'){
                        if(!openGapB){
                            difference++;
                            minLength++;
                            openGapA = false;
                            openGapB = true;
                        }
                    }
                    else if(seqA != '-' && seqB != '-'){
                        if(seqA != seqB){
                            difference++;
                        }
                        minLength++;
                        openGapA = false;
                        openGapB = false;
                    }
                }
                
                if(minLength == 0)      {    dists[h] = 1.0000;                            }
                else                    {    dists[h] = (double)difference / minLength;    }
            }
        }
        
        return dists;
    }
    catch(exception& e) {
        m->errorOut(e, "oneGapDist", "calcDist");
        exit(1);
    }
}
/***********************************************************************/

vector<int> oneGapIgnoreTermGapDist::setStarts(classifierOTU seqA, classifierOTU otu, vector<int> cols){
    try {
        vector<int> starts; starts.resize(otu.numSeqs, -1);
        
        int alignLength = cols.size();
        
        int seqAStart = 0;
        for(int i=0;i<alignLength;i++){ //for each column we want to include
            if ((seqA.otuData[cols[i]][0] != '.') && (seqA.otuData[cols[i]][0] != '-')){
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
                
                if ((thisChar == '.') || (thisChar == '-')) { } //every seq in otu is a '.' or '-' at this location, move to next column
                else { //this is a base in all locations, you are done
                    for (int k = 0; k < starts.size(); k++) {
                        if ((starts[k] == -1) && (seqA.otuData[cols[i]][0] != '.') && (seqA.otuData[cols[i]][0] != '-')) { starts[k] = i; numset++; } //any unset starts are set to this location
                    }
                    break;
                }
            }else{
                for(int j=0;j<otu.numSeqs;j++){ //for each reference
                    if((thisColumn[j] != '.') && (thisColumn[j] != '-') && (starts[j] == -1) && (seqA.otuData[cols[i]][0] != '.') && (seqA.otuData[cols[i]][0] != '-')){ //seq j hasn't set the start value and its a base
                        starts[j] = i; numset++;
                    }
                }
            }
        }
        
        return starts;
    }
    catch(exception& e) {
        m->errorOut(e, "oneGapDist", "setStarts");
        exit(1);
    }
}
/***********************************************************************/

vector<int> oneGapIgnoreTermGapDist::setEnds(classifierOTU seqA, classifierOTU otu, vector<int> cols){
    try {
        vector<int> ends; ends.resize(otu.numSeqs, -1);
        
        int alignLength = cols.size();
        
        int seqAEnd = 0;
        for(int i=alignLength-1;i>=0;i--){//for each column we want to include
            if ((seqA.otuData[cols[i]][0] != '.') && (seqA.otuData[cols[i]][0] != '-')) {
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
                
                if ((thisChar == '.') || (thisChar == '-')){ } //every seq in otu is a '.' at this location, move to next column
                else { //this is a base in all locations, you are done
                    for (int k = 0; k < ends.size(); k++) {
                        if ((ends[k] == -1) && (seqA.otuData[cols[i]][0] != '.') && (seqA.otuData[cols[i]][0] != '-')){ ends[k] = i; numset++; } //any unset starts are set to this location
                    }
                    break;
                }
            }else{
                for(int j=0;j<otu.numSeqs;j++){ //for each reference
                    if((thisColumn[j] != '.') && (thisColumn[j] != '-') && (ends[j] == -1) && (seqA.otuData[cols[i]][0] != '.') && (seqA.otuData[cols[i]][0] != '-')) { //seq j hasn't set the start value and its a base
                        ends[j] = i; numset++;
                    }
                }
            }
        }
        
        return ends;
    }
    catch(exception& e) {
        m->errorOut(e, "oneGapDist", "setEnds");
        exit(1);
    }
}

/***********************************************************************/


