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

void Rarefact::getCurve(int increment = 1, int nIters = 1000){
	try {
		RarefactionCurveData* rcd = new RarefactionCurveData();
		for(int i=0;i<displays.size();i++){
			rcd->registerDisplay(displays[i]);
		}
	
		for(int iter=0;iter<nIters;iter++){
		
			for(int i=0;i<displays.size();i++){
				displays[i]->init(label);
			}
		
			RAbundVector* lookup	= new RAbundVector(order->getNumBins());
			SAbundVector* rank	= new SAbundVector(order->getMaxRank()+1);
			random_shuffle(order->begin(), order->end());
		
			for(int i=0;i<numSeqs;i++){
			
				int binNumber = order->get(i);
				int abundance = lookup->get(binNumber);
			
				rank->set(abundance, rank->get(abundance)-1);
				abundance++;
		
				lookup->set(binNumber, abundance);
				rank->set(abundance, rank->get(abundance)+1);

				if((i == 0) || (i+1) % increment == 0){
					rcd->updateRankData(rank);
				}
			}
	
			if(numSeqs % increment != 0){
				rcd->updateRankData(rank);
			}

			for(int i=0;i<displays.size();i++){
				displays[i]->reset();
			}
		}

		for(int i=0;i<displays.size();i++){
			displays[i]->close();
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Rarefact class Function getCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Rarefact class function getCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void Rarefact::getSharedCurve(int increment = 1, int nIters = 1000){
try {
		globaldata = GlobalData::getInstance();
		SharedRarefactionCurveData* rcd = new SharedRarefactionCurveData();
		
		//register the displays
		for(int i=0;i<displays.size();i++){
			rcd->registerDisplay(displays[i]);
		}

		for(int iter=0;iter<nIters;iter++){
			//clear out old values for new iter
			lookup.clear();
		
			//create and initialize vector of sharedvectors, one for each group
			for (int i = 0; i < globaldata->Groups.size(); i++) { 
				SharedRAbundVector* temp = new SharedRAbundVector(sharedorder->getNumBins());
				temp->setLabel(sharedorder->getLabel());
				temp->setGroup(globaldata->Groups[i]);
				lookup.push_back(temp);
			}
			
			for(int i=0;i<displays.size();i++){
				displays[i]->init(label);		  
			}
		
			//sample all the members
			for(int i=0;i<numSeqs;i++){
				//get first sample
				individual chosen = sharedorder->get(i);
				int abundance; 
					
				//set info for sharedvector in chosens group
				for (int j = 0; j < lookup.size(); j++) { 
					if (chosen.group == lookup[j]->getGroup()) {
						abundance = lookup[j]->getAbundance(chosen.bin);
						lookup[j]->set(chosen.bin, (abundance + 1), chosen.group);
						break;
					}
				}
			}
			
			//randomize the groups
			random_shuffle(lookup.begin(), lookup.end());
			
			vector<SharedRAbundVector*> subset;
			//send each group one at a time
			for (int k = 0; k < lookup.size(); k++) { 
				subset.clear(); //clears out old pair of sharedrabunds
				//add in new pair of sharedrabunds
				subset.push_back(lookup[0]); subset.push_back(lookup[k]);
				
				rcd->updateSharedData(subset, k+1, numGroupComb);
				mergeVectors(lookup[0], lookup[k]);
			}

			//resets output files
			for(int i=0;i<displays.size();i++){
				displays[i]->reset();
			}
		}
		
		for(int i=0;i<displays.size();i++){
			displays[i]->close();
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Rarefact class Function getSharedCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Rarefact class function getSharedCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************************/
void Rarefact::mergeVectors(SharedRAbundVector* shared1, SharedRAbundVector* shared2) {
	try{
		for (int k = 0; k < shared1->size(); k++) {
			//merge new species into shared1
			if ((shared1->getAbundance(k) == 0) && (shared2->getAbundance(k) != 0)) {
				shared1->set(k, shared2->getAbundance(k), "combo");  //set to 'combo' since this vector now contains multiple groups
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Rarefact class Function mergeVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Rarefact class function mergeVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

