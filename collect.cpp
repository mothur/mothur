/*
 *  collect.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collect.h"


/***********************************************************************/

void Collect::getCurve(int increment = 1){
	try {
		RAbundVector* lookup = new RAbundVector(order->getNumBins());
		SAbundVector* rank	= new SAbundVector(order->getMaxRank()+1);

		CollectorsCurveData* ccd = new CollectorsCurveData();
	
		for(int i=0;i<displays.size();i++){
			ccd->registerDisplay(displays[i]); //adds a display[i] to cdd
			displays[i]->init(label);		   //sets displays label
		}									   
		for(int i=0;i<numSeqs;i++){

			int binNumber = order->get(i);
			int abundance = lookup->get(binNumber);
		
			rank->set(abundance, rank->get(abundance)-1); 
		
			abundance++;
		
			lookup->set(binNumber, abundance);
			rank->set(abundance, rank->get(abundance)+1); //increment rank(abundance)

			if((i == 0) || (i+1) % increment == 0){
				ccd->updateRankData(rank);
			}
		}
	
		if(numSeqs % increment != 0){
			ccd->updateRankData(rank);
		}
	
		for(int i=0;i<displays.size();i++){
			displays[i]->reset();
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Collect class Function getCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Collect class function getCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/***********************************************************************/
void Collect::getSharedCurve(int increment = 1){
try {
		globaldata = GlobalData::getInstance();
		vector<SharedRAbundVector*> lookup; 

		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < globaldata->gGroupmap->getNumGroups(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(sharedorder->getNumBins());
			temp->setLabel(sharedorder->getLabel());
			temp->setGroup(globaldata->gGroupmap->namesOfGroups[i]);
			temp->setGroupIndex(globaldata->gGroupmap->groupIndex[globaldata->gGroupmap->namesOfGroups[i]]);
			lookup.push_back(temp);
		}

		SharedCollectorsCurveData* ccd = new SharedCollectorsCurveData();
	
		//initialize labels for output
		//makes  'uniqueAB	 uniqueAC  uniqueBC' if your groups are A, B, C
		getGroupComb();
		groupLabel = "";
		for (int s = 0; s < groupComb.size(); s++) {
			groupLabel = groupLabel + label + groupComb[s] + "\t";
		}

		for(int i=0;i<displays.size();i++){
			ccd->registerDisplay(displays[i]); //adds a display[i] to cdd
			displays[i]->init(groupLabel);		  
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
			
			//calculate at 0 and the given increment
			if((i == 0) || (i+1) % increment == 0){
				//randomize group order
				if (globaldata->getJumble() == "1") { random_shuffle(lookup.begin(), lookup.end()); }
				//how many comparisons to make i.e. for group a, b, c = ab, ac, bc.
				int n = 1;
				for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
					for (int l = n; l < lookup.size(); l++) {
						ccd->updateSharedData(lookup[k], lookup[l], i+1, globaldata->gGroupmap->namesOfGroups.size());
					}
					n++;
				}
			}
			totalNumSeq = i+1;
		}
		
		//calculate last line if you haven't already
		if(numSeqs % increment != 0){
			//how many comparisons to make i.e. for group a, b, c = ab, ac, bc.
			int n = 1;
			for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
				for (int l = n; l < lookup.size(); l++) {
					ccd->updateSharedData(lookup[k], lookup[l], totalNumSeq, globaldata->gGroupmap->namesOfGroups.size());
				}
				n++;
			}
		}
		
		//resets output files
		for(int i=0;i<displays.size();i++){
			displays[i]->reset();
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Collect class Function getSharedCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Collect class function getSharedCurve. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************************/

void Collect::getGroupComb() {
		string group;
		
		numGroupComb = 0;
		
		int n = 1;
		for (int i = 0; i < (globaldata->gGroupmap->getNumGroups() - 1); i++) {
			for (int l = n; l < globaldata->gGroupmap->getNumGroups(); l++) {
				group = globaldata->gGroupmap->namesOfGroups[i] + globaldata->gGroupmap->namesOfGroups[l];
				groupComb.push_back(group);	
				numGroupComb++;
			}
			n++;
		}

}

/**************************************************************************************/
