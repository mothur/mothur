//
//  kimura.cpp
//  Mothur
//
//  Created by Sarah Westcott on 7/2/21.
//  Copyright © 2021 Schloss Lab. All rights reserved.
//

#include "kimura.hpp"

/***********************************************************************/

double Kimura::calcDist(Protein A, Protein B) {
    try {
        int numBases = A.getAlignLength();
        vector<AminoAcid> seqA = A.getAligned();
        vector<AminoAcid> seqB = B.getAligned();
        
        int lm = 0; int n = 0;
        
        for (int i = 0; i < numBases; i++) {
            int numA = seqA[i].getNum();
            int numB = seqB[i].getNum();
          
            if ((((long)numA <= (long)val) || ((long)numA == (long)ser)) && (((long)numB <= (long)val) || ((long)numB == (long)ser))) {
                if (numA == numB) { lm++; }
                n++;
            }
        }
        
        double p = 1 - (double)lm / n;
        double dp = 1.0 - p - 0.2 * p * p;
        
        if (dp < 0.0) {
            m->mothurOut("[WARNING]: DISTANCE BETWEEN SEQUENCES " + A.getName() + " AND " + B.getName() + " IS TOO LARGE FOR KIMURA FORMULA, setting distance to -1.0.\n");
            dist = -1.0;
        } else {
            dist = -log(dp);
        }
        
        return dist;
    }
    catch(exception& e) {
        m->errorOut(e,  "Kimura", "calcDist");
        exit(1);
    }
}
/***********************************************************************/

