//
//  f1score.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "f1score.hpp"

/***********************************************************************/
double F1Score::getValue( long long tp,  long long tn,  long long fp,  long long fn) {
    try {
        long long p = 2.0 * tp;
        long long pPrime = fn + fp;
        double f1Score = 2.0 * tp / (double) (p + pPrime);
        
        if(p + pPrime == 0)	{	f1Score = 0;	}
        
        if (isnan(f1Score) || isinf(f1Score)) { f1Score = 0; }
        
        return f1Score;
    }
    catch(exception& e) {
        m->errorOut(e, "F1Score", "getValue");
        exit(1);
    }
}
/***********************************************************************/

