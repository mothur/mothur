//
//  shannonrange.cpp
//  Mothur
//
//  Created by SarahsWork on 1/3/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "shannonrange.h"

/***********************************************************************/

EstOutput RangeShannon::getValues(SAbundVector* rank){
	try {
        data.resize(3,0);
        
        double commSize = 1e20;
        double sampleSize = rank->getNumSeqs();
        
        vector<int> freqx;
        vector<int> freqy;
        for (int i = 1; i <=rank->getMaxRank(); i++) {
            int abund = rank->get(i);
            if (abund != 0) {
                freqx.push_back(i);
                freqy.push_back(abund);
            }
        }
        
        double aux = ceil(pow((sampleSize+1), (1/(double)3)));
        double est0 = max(freqy[0]+1, aux);
        
        vector<double> ests;
        double numr = 0.0;
        double denr = 0.0;
        for (int i = 0; i < freqx.size()-1; i++) {
            
            if (m->control_pressed) { break; }
            
            if (freqx[i+1] == freqx[i]+1)   { numr = max(freqy[i+1]+1, aux);    }
            else                            { numr = aux;                       }
            
            denr = max(freqy[i], aux);
            ests.push_back((freqx[i]+1)*numr/(double)denr);
        }
        numr = aux;
        denr = max(freqy[freqy.size()-1], aux);
        ests.push_back((freqx[freqx.size()-1]+1)*numr/(double)denr);
        
        double sum = 0.0;
        for (int i = 0; i < freqy.size(); i++) {  sum += (ests[i]*freqy[i]); }
        double nfac = est0 + sum;
        est0 /= nfac;
        
        for (int i = 0; i < ests.size(); i++) {  ests[i] /= nfac;   }
        
        double abunup = 1 / commSize;
        double nbrup = est0 / abunup;
        double abunlow = ests[0];
        double nbrlow = est0 / abunlow;
        
        if (alpha == 1) {
            double sum = 0.0;
            for (int i = 0; i < freqy.size(); i++) {
                if (m->control_pressed) { break; }
                sum += (freqy[i] * ests[i] * log(ests[i]));
            }
            data[0] = -sum;
            data[1] = exp(data[0]+nbrlow*(-abunlow*log(abunlow)));
            data[2] = exp(data[0]+nbrup*(-abunup*log(abunup)));
        }else {
            for (int i = 0; i < freqy.size(); i++) {
                if (m->control_pressed) { break; }
                data[0] += (freqy[i] * (pow(ests[i],alpha)));
            }
            data[1] = pow(data[0]+nbrup*pow(abunup,alpha), (1/(1-alpha)));
            data[2] = pow(data[0]+nbrlow*pow(abunlow,alpha), (1/(1-alpha)));
        }
        
        //this calc has no data[0], just a lower and upper estimate. set data[0] to lower estimate.
        data[0] = data[1];
        if (data[1] > data[2]) { data[1] = data[2]; data[2] = data[0]; }
        data[0] = data[1];
        
       	if (isnan(data[0]) || isinf(data[0])) { data[0] = 0; }
		
		return data;
	}
	catch(exception& e) {
		m->errorOut(e, "RangeShannon", "getValues");
		exit(1);
	}
}
/***********************************************************************/
