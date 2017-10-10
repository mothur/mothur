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

int Collect::getCurve(float percentFreq = 0.01){
        try {
                RAbundVector* lookup = new RAbundVector(order->getNumBins());
                SAbundVector* rank        = new SAbundVector(order->getMaxRank()+1);

                CollectorsCurveData* ccd = new CollectorsCurveData();
        
                for(int i=0;i<displays.size();i++){
                        ccd->registerDisplay(displays[i]); //adds a display[i] to cdd
                        displays[i]->init(label);                   //sets displays label
                }   
				
				//convert freq percentage to number
				int increment = 1;
				if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
				else { increment = percentFreq;  }
																						                                                                        
                for(int i=0;i<numSeqs;i++){
						
						if (m->getControl_pressed()) { delete lookup; delete rank; delete ccd;  return 1;  }
						
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
				
				return 0;
        }
        catch(exception& e) {
			m->errorOut(e, "Collect", "getCurve");
			exit(1);
        }
}

/***********************************************************************/
int Collect::getSharedCurve(float percentFreq = 0.01){
try {
                vector<SharedRAbundVector*> lookup;
                map<string, int> indexLookup;
				vector<SharedRAbundVector*> subset;

                //create and initialize vector of sharedvectors, one for each group
				vector<string> mGroups = m->getGroups();
                for (int i = 0; i < mGroups.size(); i++) {
                        SharedRAbundVector* temp = new SharedRAbundVector(sharedorder->getNumBins());
                        temp->setLabel(sharedorder->getLabel());
                        temp->setGroup(mGroups[i]);
                        indexLookup[mGroups[i]] = i;
						lookup.push_back(temp);
                }
	
                SharedCollectorsCurveData* ccd = new SharedCollectorsCurveData();
        
                //initialize labels for output
                //makes  'uniqueAB         uniqueAC  uniqueBC' if your groups are A, B, C
                getGroupComb();
				
                for(int i=0;i<displays.size();i++){
                        ccd->registerDisplay(displays[i]); //adds a display[i] to cdd
						bool hasLciHci = displays[i]->hasLciHci();
						groupLabel = "";
						for (int s = 0; s < groupComb.size(); s++) {
							if (hasLciHci) {  groupLabel = groupLabel + label + groupComb[s] + "\t" + label + groupComb[s] + "lci\t" + label + groupComb[s] + "hci\t"; }
							else{  groupLabel = groupLabel + label + groupComb[s] + "\t";  }
						}

						string groupLabelAll = groupLabel + label + "all\t"; 
						if ((displays[i]->isCalcMultiple() ) && (displays[i]->getAll() )) {   displays[i]->init(groupLabelAll); }
						else {  displays[i]->init(groupLabel);  }           
                }
                
				//convert freq percentage to number
				int increment = 1;
				if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
				else { increment = percentFreq;  }
				
                //sample all the members
                for(int i=0;i<numSeqs;i++){
				
						if (m->getControl_pressed()) { for (int j = 0; j < lookup.size(); j++) {  delete lookup[j]; } delete ccd;  return 1;  }
						
                        //get first sample
                        individual chosen = sharedorder->get(i);
                        int abundance = lookup[indexLookup[chosen.group]]->get(chosen.bin);
                        lookup[indexLookup[chosen.group]]->set(chosen.bin, (abundance + 1));
                    
	
                        //calculate at 0 and the given increment
                        if((i == 0) || (i+1) % increment == 0){

								//how many comparisons to make i.e. for group a, b, c = ab, ac, bc.

                                int n = 1;
								bool pair = true;
                                for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
                                        for (int l = n; l < lookup.size(); l++) {
												subset.clear(); //clear out old pair of sharedrabunds
												//add new pair of sharedrabund vectors
												subset.push_back(lookup[k]); subset.push_back(lookup[l]);
											
												//load subset with rest of lookup for those calcs that need everyone to calc for a pair
												for (int w = 0; w < lookup.size(); w++) {
													if ((w != k) && (w != l)) { subset.push_back(lookup[w]); }
												}
						
                                                ccd->updateSharedData(subset, i+1, m->getNumGroups(), pair);
                                        }
                                        n++;
                                }
							
								//if this is a calculator that can do multiples then do them
								pair = false;
								ccd->updateSharedData(lookup, i+1, m->getNumGroups(), pair); 
							
                        }
                        totalNumSeq = i+1;
                }
                
                //calculate last label if you haven't already
                if(numSeqs % increment != 0){
                        //how many comparisons to make i.e. for group a, b, c = ab, ac, bc.
                        int n = 1;
						bool pair = true;
                        for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
                                for (int l = n; l < lookup.size(); l++) {
										subset.clear(); //clear out old pair of sharedrabunds
										//add new pair of sharedrabund vectors
										subset.push_back(lookup[k]); subset.push_back(lookup[l]);
									
										//load subset with rest of lookup for those calcs that need everyone to calc for a pair
										for (int w = 0; w < lookup.size(); w++) {
											if ((w != k) && (w != l)) { subset.push_back(lookup[w]); }
										}
									
										ccd->updateSharedData(subset, totalNumSeq, m->getNumGroups(), pair);
                                }
                                n++;
                        }
						//if this is a calculator that can do multiples then do them
						pair = false;
						ccd->updateSharedData(lookup, totalNumSeq, m->getNumGroups(), pair); 
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
				
				return 0;

        }
        catch(exception& e) {
                m->errorOut(e, "Collect", "getSharedCurve");
				exit(1);
        }
}

/**************************************************************************************/

void Collect::getGroupComb() {
	string group;
                
	numGroupComb = 0;
                
	int n = 1;
	vector<string> mGroups = m->getGroups();
	for (int i = 0; i < (m->getNumGroups() - 1); i++) {
		for (int l = n; l < m->getNumGroups(); l++) {
			group = mGroups[i] + mGroups[l];
			groupComb.push_back(group);        
			numGroupComb++;
		}
		n++;
	}

}

/**************************************************************************************/
