/*
 *  rarefact.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefact.h"


/**************************************************************************************************/
struct singleRarefactData {
    
    long long nIters, numSeqs;
    MothurOut* m;
    Utils util;
    OrderVector order;
    set<int> ends;
    vector<Display*> displays;
    string label;
    int increment;
    
    singleRarefactData(){}
    singleRarefactData(long long st, Utils u, OrderVector o, set<int> ed, vector<Display*>& dis, string l, long long ns, int inc) {
        m = MothurOut::getInstance();
        nIters = st;
        util = u;
        order = o;
        ends = ed;
        displays = dis;
        label = l;
        numSeqs = ns;
        increment = inc;
    }
};
/***********************************************************************/
int singleDriver(singleRarefactData* params){
    try {
        
        RarefactionCurveData rcd; rcd.registerDisplays(params->displays);
        
        for(int iter=0;iter<params->nIters;iter++){
            
            for(int i=0;i<params->displays.size();i++){ params->displays[i]->init(params->label); }
            
            RAbundVector lookup(params->order.getNumBins());
            SAbundVector* rank	= new SAbundVector(params->order.getMaxRank()+1);
            params->util.mothurRandomShuffle(params->order);
            
            for(int i=0;i<params->numSeqs;i++){
                
                if (params->m->getControl_pressed()) { delete rank; return 0;  }
                
                int binNumber = params->order.get(i);
                int abundance = lookup.get(binNumber);
                
                rank->set(abundance, rank->get(abundance)-1);
                abundance++;
                
                lookup.set(binNumber, abundance);
                rank->set(abundance, rank->get(abundance)+1);
                
                if((i == 0) || ((i+1) % params->increment == 0) || (params->ends.count(i+1) != 0)){ rcd.updateRankData(rank); }
            }
            
            if((params->numSeqs % params->increment != 0) || (params->ends.count(params->numSeqs) != 0)){ rcd.updateRankData(rank); }
            
            for(int i=0;i<params->displays.size();i++){ params->displays[i]->reset(); }
            
            delete rank;
        }
        
        return 0;
    }
    catch(exception& e) {
        params->m->errorOut(e, "Rarefact", "singleDriver");
        exit(1);
    }
}

/***********************************************************************/

int Rarefact::getCurve(float percentFreq = 0.01, int nIters = 1000){
	try {
		//convert freq percentage to number
		int increment = 1;
		if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
		else { increment = percentFreq;  }
        
        vector<int> lines;
        if (processors > (nIters)) { processors = nIters; }
        
        //figure out how many sequences you have to process
        int numItersPerProcessor = nIters / processors;
        for (int i = 0; i < processors; i++) {
            if(i == (processors - 1)){	numItersPerProcessor = (nIters) - i * numItersPerProcessor; 	}
            lines.push_back(numItersPerProcessor);
        }
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<singleRarefactData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            
            //make copy of order so we don't get access violations
            OrderVector newOrder(order);
            singleRarefactData* dataBundle = new singleRarefactData(lines[i+1], util, newOrder, ends, displays, label, numSeqs, increment);
            
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(singleDriver, dataBundle));
        }
        
        //make copy of lookup so we don't get access violations
        OrderVector newOrder(order);
        singleRarefactData* dataBundle = new singleRarefactData(lines[0], util, newOrder, ends, displays, label, numSeqs, increment);
        singleDriver(dataBundle);
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            
            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;

		for(int i=0;i<displays.size();i++){ displays[i]->close(); }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Rarefact", "getCurve");
		exit(1);
	}
}

/***********************************************************************/
int Rarefact::getSharedCurve(float percentFreq = 0.01, int nIters = 1000){
try {
		SharedRarefactionCurveData rcd;
		
		label = lookup[0]->getLabel();
		
		//register the displays
		rcd.registerDisplays(displays);
		
		//if jumble is false all iters will be the same
		if (!jumble)  {  nIters = 1;  }
		
		//convert freq percentage to number
		int increment = 1;
		if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
		else { increment = percentFreq;  }
		
		for(int iter=0;iter<nIters;iter++){
		
			for(int i=0;i<displays.size();i++){ displays[i]->init(label);	 }
			
            //randomize the groups
			if (jumble)  { util.mothurRandomShuffle(lookup); }
			
			//make merge the size of lookup[0]
			SharedRAbundVector* merge = new SharedRAbundVector(lookup[0]->getNumBins());
			
			//make copy of lookup zero
			for(int i = 0; i<lookup[0]->getNumBins(); i++) {  merge->set(i, lookup[0]->get(i)); }
			
			vector<SharedRAbundVector*> subset;
			//send each group one at a time
			for (int k = 0; k < lookup.size(); k++) { 
				if (m->getControl_pressed()) {  delete merge;  return 0;  }
				
				subset.clear(); //clears out old pair of sharedrabunds
				//add in new pair of sharedrabunds
				subset.push_back(merge); subset.push_back(lookup[k]);
				
				rcd.updateSharedData(subset, k+1, numGroupComb); //
				mergeVectors(merge, lookup[k]);
			}

			//resets output files
			for(int i=0;i<displays.size();i++){ displays[i]->reset(); }
			
			delete merge;
		}
		
		for(int i=0;i<displays.size();i++){ displays[i]->close(); }
		
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

