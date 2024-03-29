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
                
               // printf("l = %ld, nb1 = %ld, nb2 = %ld\n", l, nb1[i], nb2[i]);
               // printf("l = %ld, eig[m] = %f, prob[m][nb1 - 1] = %f, prob[m][nb2 - 1] = %f\n", l, eigs[l], probs[l][nb1[i] - 1], probs[l][nb2[i] - 1]);

                
                q = probs[l][nb1[i]-1] * probs[l][nb2[i]-1] * elambdat;
                p += q;
                
                dp += eigs[l] * q;
                
                double TEMP = eigs[l];
                
                d2p += TEMP * TEMP * q;
            }
        }
        
        //printf("p = %f, q = %f, tt = %f\n", p, q, tt);
       // printf("dp = %f, d2p = %f\n", dp, d2p);
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
double DistCalc::makeDists(Protein A, Protein B, double eigs[20], double probs[20][20]){
    try {
        int numBases = A.getAlignLength();
        vector<AminoAcid> seqA = A.getAligned();
        vector<AminoAcid> seqB = B.getAligned();
        
        bool inf = false; bool neginfinity = false; bool overlap = false;
        double delta, lnlike, slope, curv, tt;
        tt = 0.1; delta = tt / 2.0;
        
        for (int l = 0; l < 20; l++) {
            
            //reset for this attempt
            lnlike = 0.0; slope = 0.0; curv = 0.0; neginfinity = false; overlap = false;
            double oldweight = 1.0;
            
            if (m->getControl_pressed()) { break; }
            
            for (int i = 0; i < numBases; i++) {
                int numA = seqA[i].getNum();
                int numB = seqB[i].getNum();
                
                if (numA != stop && numA != del && numA != quest && numA != unk &&
                    numB != stop && numB != del && numB != quest && numB != unk) {
            
                    double p = 0.0; double dp = 0.0; double d2p = 0.0; overlap = true;
                    
                    vector<int> numAs; vector<int> numBs;
                    if (numA != asx && numA != glx && numB != asx && numB != glx) {
                        if (numA < ser2) { numA++; }
                        if (numB < ser2) { numB++; }
                        numAs.push_back(numA); numBs.push_back(numB); //+1 avoid 0
                    }else {
                        fillNums(numAs, numBs, numA, numB);
                    }
                    
                    predict(numAs, numBs, p, dp, d2p, tt, eigs, probs);
                    
                    if (p <= 0.0) {
                        neginfinity = true;
                    }else {
                        slope += oldweight*dp / p;
                        curv += oldweight*(d2p / p - dp * dp / (p * p));
                        
                        //printf("%ld:%ld, dp = %f, p = %f, d2p = %f\n", l, i, dp, p, d2p);

                        //printf("%ld:%ld, slope = %f, curv = %f, oldweight[i] = %ld\n", l, i, slope, curv, oldweight);
                    }
                }//endif stop
            }//endif bases
            
            if (!overlap) {
                tt = -1.0; l += 20; inf = true;
            }else if (!neginfinity) {
                if (curv < 0.0) {
                    tt -= slope / curv;
                    if (tt > 10000.0) { tt = -1.0; l += 20; inf = true;  }
                }else {
                    if ((slope > 0.0 && delta < 0.0) || (slope < 0.0 && delta > 0.0)) { delta /= -2; }
                    tt += delta;
                }
            }else {
                delta /= -2;
                tt += delta;
            }
            
            if (tt < 0.00001 && !inf) { tt = 0.00001; }
            
        }//endif attempts
        
        dist = tt;
        //exit(1);
        
        return dist;
    }
    catch(exception& e) {
        m->errorOut(e,  "DistCalc", "makeDists");
        exit(1);
    }
}
/***********************************************************************/
