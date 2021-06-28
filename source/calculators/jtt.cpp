//
//  jtt.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/26/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "jtt.hpp"

/***********************************************************************/
double JTT::calcDist(Protein A, Protein B){
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
                    
                    predict(numAs, numBs, p, dp, d2p, tt);
                    
                    if (p <= 0.0) {
                        neginfinity = true;
                    }else {
                        slope += oldweight*dp / p;
                        curv += oldweight*(d2p / p - dp * dp / (p * p));
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
        
        return dist;
    }
    catch(exception& e) {
        m->errorOut(e,  "JTT", "calcDist");
        exit(1);
    }
}
/***********************************************************************/
//nb1 and nb2 have size 1, unless amino acid = B or Z
void JTT::predict(vector<int> nb1, vector<int> nb2, double& p, double& dp, double& d2p, double& tt){
    try {
        double q;
        
        for (int i = 0; i < nb1.size(); i++) {
            
            for (int l = 0; l < 20; l++) {
              
                double elambdat = exp(tt * jtteigs[l]);
                
                q = jttprobs[l][nb1[i]-1] * jttprobs[l][nb2[i]-1] * elambdat;
                p += q;
                
                dp += jtteigs[l] * q;
                
                double TEMP = jtteigs[l];
                
                d2p += TEMP * TEMP * q;
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e,  "JTT", "predict");
        exit(1);
    }
}
/***********************************************************************/
//nb1 and nb2 have size 1, unless amino acid = B or Z
void JTT::fillNums(vector<int>& numAs, vector<int>& numBs, int numA, int numB){
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
        m->errorOut(e,  "JTT", "fillNums");
        exit(1);
    }
}
/***********************************************************************/
