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
                SAbundVector* rank        = new SAbundVector(order->getMaxRank()+1);

                CollectorsCurveData* ccd = new CollectorsCurveData();
        
                for(int i=0;i<displays.size();i++){
                        ccd->registerDisplay(displays[i]); //adds a display[i] to cdd
                        displays[i]->init(label);                   //sets displays label
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
				
				delete lookup;
				delete rank;
				delete ccd;
        }
        catch(exception& e) {
			errorOut(e, "Collect", "getCurve");
			exit(1);
        }
}

/***********************************************************************/
void Collect::getSharedCurve(int increment = 1){
try {
                globaldata = GlobalData::getInstance();
                vector<SharedRAbundVector*> lookup; 
				vector<SharedRAbundVector*> subset;

                //create and initialize vector of sharedvectors, one for each group
                for (int i = 0; i < globaldata->Groups.size(); i++) { 
                        SharedRAbundVector* temp = new SharedRAbundVector(sharedorder->getNumBins());
                        temp->setLabel(sharedorder->getLabel());
                        temp->setGroup(globaldata->Groups[i]);
                        temp->setGroupIndex(globaldata->gGroupmap->groupIndex[globaldata->Groups[i]]);
                        lookup.push_back(temp);
                }

                SharedCollectorsCurveData* ccd = new SharedCollectorsCurveData();
        
                //initialize labels for output
                //makes  'uniqueAB         uniqueAC  uniqueBC' if your groups are A, B, C
                getGroupComb();
                groupLabel = ""; 
                for (int s = 0; s < groupComb.size(); s++) {
                        groupLabel = groupLabel + label + groupComb[s] + "\t";
                }
				
				//for multi displays
				string groupLabelAll = groupLabel + label + "all\t";
				
                for(int i=0;i<displays.size();i++){
                        ccd->registerDisplay(displays[i]); //adds a display[i] to cdd
						if (displays[i]->isCalcMultiple() == true)  {   displays[i]->init(groupLabelAll); }
						else {  displays[i]->init(groupLabel);  }           
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

								//how many comparisons to make i.e. for group a, b, c = ab, ac, bc.

                                int n = 1;
                                for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
                                        for (int l = n; l < lookup.size(); l++) {
												subset.clear(); //clear out old pair of sharedrabunds
												//add new pair of sharedrabund vectors
												subset.push_back(lookup[k]); subset.push_back(lookup[l]);
                                                ccd->updateSharedData(subset, i+1, globaldata->Groups.size());
                                        }
                                        n++;
                                }
								//if this is a calculator that can do multiples then do them
								ccd->updateSharedData(lookup, i+1, globaldata->Groups.size()); 
                        }
                        totalNumSeq = i+1;
                }
                
                //calculate last line if you haven't already
                if(numSeqs % increment != 0){
                        //how many comparisons to make i.e. for group a, b, c = ab, ac, bc.
                        int n = 1;
                        for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
                                for (int l = n; l < lookup.size(); l++) {
										subset.clear(); //clear out old pair of sharedrabunds
										//add new pair of sharedrabund vectors
										subset.push_back(lookup[k]); subset.push_back(lookup[l]);
										ccd->updateSharedData(subset, totalNumSeq, globaldata->Groups.size());
                                }
                                n++;
                        }
						//if this is a calculator that can do multiples then do them
						ccd->updateSharedData(lookup, totalNumSeq, globaldata->Groups.size()); 
                }
                
                //resets output files
                for(int i=0;i<displays.size();i++){
                        displays[i]->reset();
                }
				
				//memory cleanup
				delete ccd;
				for (int i = 0; i < lookup.size(); i++) {
					delete lookup[i];
				}

        }
        catch(exception& e) {
                errorOut(e, "Collect", "getSharedCurve");
				exit(1);
        }
}

/**************************************************************************************/

void Collect::getGroupComb() {
	string group;
                
	numGroupComb = 0;
                
	int n = 1;
	for (int i = 0; i < (globaldata->Groups.size() - 1); i++) {
		for (int l = n; l < globaldata->Groups.size(); l++) {
			group = globaldata->Groups[i] + globaldata->Groups[l];
			groupComb.push_back(group);        
			numGroupComb++;
		}
		n++;
	}

}

/**************************************************************************************/
