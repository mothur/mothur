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
//nb1 and nb2 have size 1, unless amino acid = B or Z
void DistCalc::predict(vector<int> nb1, vector<int> nb2, double& p, double& dp, double& d2p, double& tt, double eigs[20], double probs[20][20]){
    try {
        double q;
        
        for (int i = 0; i < nb1.size(); i++) {
            
            for (int l = 0; l < 20; l++) {
              
                double elambdat = exp(tt * eigs[l]);
                
                q = probs[l][nb1[i]-1] * probs[l][nb2[i]-1] * elambdat;
                p += q;
                
                dp += eigs[l] * q;
                
                double TEMP = eigs[l];
                
                d2p += TEMP * TEMP * q;
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e,  "DistCalc", "predict");
        exit(1);
    }
}
/***********************************************************************/
//nb1 and nb2 have size 1, unless amino acid = B or Z
void DistCalc::fillNums(vector<int>& numAs, vector<int>& numBs, int numA, int numB){
    try {
        
        if (numA == asx) { //B asn or asp (3 or 4)
          if (numB == asx) { //B asn or asp
              numAs.push_back(3); numBs.push_back(3); //asn, asn
              numAs.push_back(3); numBs.push_back(4); //asn, asp
              numAs.push_back(4); numBs.push_back(3); //asp, asn
              numAs.push_back(4); numBs.push_back(4); //asp, asp
          }else {
              if (numB == glx) { //Z gln or glu (6 or 7)
                  numAs.push_back(3); numBs.push_back(6); //asn, gln
                  numAs.push_back(3); numBs.push_back(7); //asn, glu
                  numAs.push_back(4); numBs.push_back(6); //asp, gln
                  numAs.push_back(4); numBs.push_back(7); //asp, glu
              }else {
                  if (numB < ser2) { numB++; }
                  numAs.push_back(3); numBs.push_back(numB); //asn, numB
                  numAs.push_back(4); numBs.push_back(numB); //asp, numB
              }
          }
        }else {
            if (numA == glx) { //Z gln or glu (6 or 7)
                if (numB == asx) { //B asn or asp
                    numAs.push_back(6); numBs.push_back(3); //gln, asn
                    numAs.push_back(6); numBs.push_back(4); //gln, asp
                    numAs.push_back(7); numBs.push_back(3); //glu, asn
                    numAs.push_back(7); numBs.push_back(4); //glu, asp
                }else {
                    if (numB == glx) { //Z gln or glu (6 or 7)
                        numAs.push_back(6); numBs.push_back(6); //gln, gln
                        numAs.push_back(6); numBs.push_back(7); //gln, glu
                        numAs.push_back(7); numBs.push_back(6); //glu, gln
                        numAs.push_back(7); numBs.push_back(7); //glu, glu
                    }else {
                        if (numB < ser2) { numB++; }
                        numAs.push_back(6); numBs.push_back(numB); //gln, numB
                        numAs.push_back(7); numBs.push_back(numB); //glu, numB
                    }
                }
            }else {
                if (numA < ser2) { numA++; }
                if (numB == asx) { //B asn or asp
                    numAs.push_back(numA); numBs.push_back(3); //numA, asn
                    numAs.push_back(numA); numBs.push_back(4); //numA, asp
                    numAs.push_back(numA); numBs.push_back(3); //numA, asn
                    numAs.push_back(numA); numBs.push_back(4); //numA, asp
                }else if (numB == glx) { //Z gln or glu (6 or 7)
                    numAs.push_back(numA); numBs.push_back(6); //numA, gln
                    numAs.push_back(numA); numBs.push_back(7); //numA, glu
                    numAs.push_back(numA); numBs.push_back(6); //numA, gln
                    numAs.push_back(numA); numBs.push_back(7); //numA, glu
                }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e,  "DistCalc", "fillNums");
        exit(1);
    }
}
/***********************************************************************/
