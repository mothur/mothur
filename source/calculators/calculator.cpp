//
//  calculator.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/21/20.
//  Copyright © 2020 Schloss Lab. All rights reserved.
//

#include "calculator.h"

/***********************************************************************/
int DistCalc::setStart(string seqA, string seqB) {
    try {
        int start = 0;
        int alignLength = seqA.length();
        
        for(int i=0;i<alignLength;i++){
            if((seqA[i] != '.' || seqB[i] != '.')){ //one of you is not a terminal gap
                start = i;
                break;
            }
        }
        
        return start;
    }
    catch(exception& e) {
        m->errorOut(e, "DistCalc", "setStart");
        exit(1);
    }
}
/***********************************************************************/
int DistCalc::setEnd(string seqA, string seqB) {
    try {
        int end = 0;
        int alignLength = seqA.length();
        
        for(int i=alignLength-1;i>=0;i--){
            if((seqA[i] != '.' || seqB[i] != '.')){ //one of you is not a terminal gap
                end = i;
                break;
            }
        }
        return end;
    }
    catch(exception& e) {
        m->errorOut(e, "DistCalc", "setEnd");
        exit(1);
    }
}
/***********************************************************************/
// this assumes that sequences start and end with '.'s instead of'-'s.
int DistCalc::setStartIgnoreTermGap(string seqA, string seqB, bool& overlap) {
    try {
        
        int start = 0;
        int alignLength = seqA.length();
        
        for(int i=0;i<alignLength;i++){
            if(seqA[i] != '.' && seqB[i] != '.' && seqA[i] != '-' && seqB[i] != '-' ){ //skip leading gaps
                start = i;
                
                overlap = true;
                break;
            }
        }
        
        return start;
    }
    catch(exception& e) {
        m->errorOut(e, "DistCalc", "setStartIgnoreTermGap");
        exit(1);
    }
}
/***********************************************************************/
// this assumes that sequences start and end with '.'s instead of'-'s.
int DistCalc::setEndIgnoreTermGap(string seqA, string seqB, bool& overlap) {
    try {
        
        int end = 0;
        int alignLength = seqA.length();
        
        for(int i=alignLength-1;i>=0;i--){
            if(seqA[i] != '.' && seqB[i] != '.' && seqA[i] != '-' && seqB[i] != '-' ){ //ignore terminal gaps
                end = i;
                
                overlap = true;
                break;
            }
        }
        
        return end;
    }
    catch(exception& e) {
        m->errorOut(e, "DistCalc", "setEndIgnoreTermGap");
        exit(1);
    }
}
/***********************************************************************/

vector<int> DistCalc::setStartsIgnoreTermGap(classifierOTU seqA, classifierOTU otu, vector<int> cols){
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
        m->errorOut(e, "DistCalc", "setStartsIgnoreTermGap");
        exit(1);
    }
}
/***********************************************************************/

vector<int> DistCalc::setEndsIgnoreTermGap(classifierOTU seqA, classifierOTU otu, vector<int> cols){
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
        m->errorOut(e, "DistCalc", "setEndsIgnoreTermGap");
        exit(1);
    }
}
/***********************************************************************/

vector<int> DistCalc::setStarts(classifierOTU seqA, classifierOTU otu, vector<int> cols){
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
        m->errorOut(e, "DistCalc", "setStarts");
        exit(1);
    }
}
/***********************************************************************/

vector<int> DistCalc::setEnds(classifierOTU seqA, classifierOTU otu, vector<int> cols){
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
        m->errorOut(e, "DistCalc", "setEnds");
        exit(1);
    }
}

/***********************************************************************/
