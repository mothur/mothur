//
//  pmb.cpp
//  Mothur
//
//  Created by Sarah Westcott on 6/29/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "pmb.hpp"

/***********************************************************************/
double PMB::calcDist(Protein A, Protein B){
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
                    
                    predict(numAs, numBs, p, dp, d2p, tt, pmbeigs, pmbprobs);
                    
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
