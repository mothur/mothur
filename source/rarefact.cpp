/*
 *  rarefact.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefact.h"
//#include "ordervector.hpp"

/***********************************************************************/

int Rarefact::getCurve(float percentFreq = 0.01, int nIters = 1000){
	try {
		RarefactionCurveData* rcd = new RarefactionCurveData();
		for(int i=0;i<displays.size();i++){ rcd->registerDisplay(displays[i]); }
		
		//convert freq percentage to number
		int increment = 1;
		if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
		else { increment = percentFreq;  }	
		
        driver(rcd, increment, nIters);

		for(int i=0;i<displays.size();i++){ displays[i]->close(); }
		
		delete rcd;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Rarefact", "getCurve");
		exit(1);
	}
}
/***********************************************************************/
int Rarefact::driver(RarefactionCurveData* rcd, int increment, int nIters = 1000){
	try {
			
		for(int iter=0;iter<nIters;iter++){
		
			for(int i=0;i<displays.size();i++){ displays[i]->init(label); }
		
			RAbundVector* lookup	= new RAbundVector(order.getNumBins());
			SAbundVector* rank	= new SAbundVector(order.getMaxRank()+1);
			util.mothurRandomShuffle(order);
		
			for(int i=0;i<numSeqs;i++){
			
				if (m->getControl_pressed()) { delete lookup; delete rank; delete rcd; return 0;  }
			
				int binNumber = order.get(i);
				int abundance = lookup->get(binNumber);
			
				rank->set(abundance, rank->get(abundance)-1);
				abundance++;
		
				lookup->set(binNumber, abundance);
				rank->set(abundance, rank->get(abundance)+1);

				if((i == 0) || ((i+1) % increment == 0) || (ends.count(i+1) != 0)){ rcd->updateRankData(rank); }
			}
	
			if((numSeqs % increment != 0) || (ends.count(numSeqs) != 0)){ rcd->updateRankData(rank); }

			for(int i=0;i<displays.size();i++){ displays[i]->reset(); }
			
			delete lookup;
			delete rank;
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Rarefact", "driverRarefact");
		exit(1);
	}
}
/***********************************************************************/
int Rarefact::getSharedCurve(float percentFreq = 0.01, int nIters = 1000){
try {
		SharedRarefactionCurveData* rcd = new SharedRarefactionCurveData();
		
		label = lookup[0]->getLabel();
		
		//register the displays
		for(int i=0;i<displays.size();i++){
			rcd->registerDisplay(displays[i]);
		}
		
		//if jumble is false all iters will be the same
		if (!jumble)  {  nIters = 1;  }
		
		//convert freq percentage to number
		int increment = 1;
		if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
		else { increment = percentFreq;  }
		
		for(int iter=0;iter<nIters;iter++){
		
			for(int i=0;i<displays.size();i++){
				displays[i]->init(label);		  
			}
			
            //randomize the groups
			if (jumble)  { util.mothurRandomShuffle(lookup); }
			
			//make merge the size of lookup[0]
			SharedRAbundVector* merge = new SharedRAbundVector(lookup[0]->getNumBins());
			
			//make copy of lookup zero
			for(int i = 0; i<lookup[0]->getNumBins(); i++) {  merge->set(i, lookup[0]->get(i)); }
			
			vector<SharedRAbundVector*> subset;
			//send each group one at a time
			for (int k = 0; k < lookup.size(); k++) { 
				if (m->getControl_pressed()) {  delete merge; delete rcd; return 0;  }
				
				subset.clear(); //clears out old pair of sharedrabunds
				//add in new pair of sharedrabunds
				subset.push_back(merge); subset.push_back(lookup[k]);
				
				rcd->updateSharedData(subset, k+1, numGroupComb);
				mergeVectors(merge, lookup[k]);
			}

			//resets output files
			for(int i=0;i<displays.size();i++){
				displays[i]->reset();
			}
			
			delete merge;
		}
		
		for(int i=0;i<displays.size();i++){
			displays[i]->close();
		}
		
		delete rcd;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Rarefact", "getSharedCurve");
		exit(1);
	}
}

/**************************************************************************************/
void Rarefact::mergeVectors(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
	try{
		for (int k = 0; k < shared1->getNumBins(); k++) {
			//merge new species into shared1
			shared1->set(k, (shared1->get(k) + shared2->get(k)));  //set to 'combo' since this vector now contains multiple groups
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Rarefact", "mergeVectors");
		exit(1);
	}
}

