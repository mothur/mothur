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
        
        int alignLength = cols.size();
        
        if (setCutoff) {
            
            vector<int> ends = setEnds(seq, otu, cols);
            
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
            
        }else{
            for (int h = 0; h < otu.numSeqs; h++) {
                
                if (m->getControl_pressed()) { break; }
                
                int minLength = 0;
                int difference = 0;
                bool openGapA = false;
                bool openGapB = false;
                
                for(int i=starts[h];i<alignLength;i++){
                    
                    char seqA = seq.otuData[cols[i]][0];
                    vector<char> otuChars = otu.otuData[cols[i]];
                    char seqB = otuChars[0]; //assume column if identical
                    if (otuChars.size() == otu.numSeqs) {  seqB = otuChars[h]; }
                    
                    
                    if((seqA == '-' && seqB == '-') || (seqA == '.' && seqB == '-') || (seqA == '-' && seqB == '.')){    ;    }
                    else if(seqA == '.' && seqB == '.'){ break; }
                    
                    else if(seqB != '-' && (seqA == '-' || seqA == '.')){
                        if(!openGapA){
                            difference++;
                            openGapA = true;
                            openGapB = false;
                            minLength++;
                        }
                    }
                    else if(seqA != '-' && (seqB == '-' || seqB == '.')){
                        if(!openGapB){
                            difference++;
                            openGapA = false;
                            openGapB = true;
                            minLength++;
                        }
                    }
                    else if(seqA != '-' && seqB != '-'){
                        minLength++;
                        openGapA = false;
                        openGapB = false;
                        if(seqA != seqB){ difference++; }
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

vector<int> oneGapDist::setStarts(classifierOTU seqA, classifierOTU otu, vector<int> cols){
    try {
        vector<int> starts; starts.resize(otu.numSeqs, 0);
        
        int alignLength = cols.size();
        
        int seqAStart = 0;
        for(int i=0;i<alignLength;i++){ //for each column we want to include
            if (seqA.otuData[cols[i]][0] != '.') {
                seqAStart = i; break;
            }
        }
        
        //set start positions
        int numset = 0;
        for(int i=0;i<alignLength;i++){ //for each column we want to include
            
            if (seqAStart <= i) { //our query seq starts before this point so set rest of unset starts to query start
                for (int k = 0; k < starts.size(); k++) {
                    if (starts[k] == 0) { starts[k] = seqAStart; numset++; }
                }
                break;
            }else if(numset == otu.numSeqs) { break; }
            
            vector<char> thisColumn = otu.otuData[cols[i]];
            if (thisColumn.size() != otu.numSeqs) { //all seqs at this spot are identical
                
                char thisChar = thisColumn[0];
                
                if (thisChar == '.') { } //every seq in otu is a '.' at this location, move to next column
                else { //this is a base in all locations, you are done
                    for (int k = 0; k < starts.size(); k++) {
                        if (starts[k] == 0) { starts[k] = i; numset++; } //any unset starts are set to this location
                    }
                    break;
                }
            }else{
                for(int j=0;j<otu.numSeqs;j++){ //for each reference
                    if((thisColumn[j] != '.') && (starts[j] == 0)){ //seq j hasn't set the start value and its a base
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

vector<int> oneGapDist::setEnds(classifierOTU seqA, classifierOTU otu, vector<int> cols){
    try {
        vector<int> ends; ends.resize(otu.numSeqs, 0);
        
        int alignLength = cols.size();
        
        int seqAEnd = 0;
        for(int i=alignLength-1;i>=0;i--){//for each column we want to include
            if (seqA.otuData[cols[i]][0] != '.') {
                seqAEnd = i; break;
            }
        }
        
        //set start positions
        int numset = 0;
        for(int i=alignLength-1;i>=0;i--){ //for each column we want to include
            
            if (seqAEnd <= i) { //our query seq starts before this point so set rest of unset starts to query start
                for (int k = 0; k < ends.size(); k++) {
                    if (ends[k] == 0) { ends[k] = seqAEnd; numset++; }
                }
                break;
            }else if(numset == otu.numSeqs) { break; }
            
            vector<char> thisColumn = otu.otuData[cols[i]];
            if (thisColumn.size() != otu.numSeqs) { //all seqs at this spot are identical
                
                char thisChar = thisColumn[0];
                
                if (thisChar == '.') { } //every seq in otu is a '.' at this location, move to next column
                else { //this is a base in all locations, you are done
                    for (int k = 0; k < ends.size(); k++) {
                        if (ends[k] == 0) { ends[k] = i; numset++; } //any unset starts are set to this location
                    }
                    break;
                }
            }else{
                for(int j=0;j<otu.numSeqs;j++){ //for each reference
                    if((thisColumn[j] != '.') && (ends[j] == 0)){ //seq j hasn't set the start value and its a base
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



