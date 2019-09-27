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
            RAbundVector rabund(order->getNumBins());
            SAbundVector rank(order->getMaxRank()+1);

            //sets displays label
            for(int i=0;i<displays.size();i++){ displays[i]->init(label); }
            
            CollectorsCurveData ccd; ccd.registerDisplays(displays);
            
            //convert freq percentage to number
            int increment = 1;
            if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
            else { increment = percentFreq;  }
            
            for(int i=0;i<numSeqs;i++){
                
                if (m->getControl_pressed()) {   return 1;  }
                
                int binNumber = order->get(i);
                int abundance = rabund.get(binNumber);
                
                rank.set(abundance, rank.get(abundance)-1);
                
                abundance++;
                
                rabund.set(binNumber, abundance);
                rank.set(abundance, rank.get(abundance)+1); //increment rank(abundance)
                
                if((i == 0) || (i+1) % increment == 0){ ccd.updateRankData(rank); }
            }
            
            if(numSeqs % increment != 0){ ccd.updateRankData(rank); }
            
            for(int i=0;i<displays.size();i++){ displays[i]->reset(); }
            
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
        
        //create and initialize vector of sharedvectors, one for each group
        vector<string> groups = sharedorder->getGroups();
        int numGroups = groups.size();
        for (int i = 0; i < groups.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(sharedorder->getMaxRank()+1);
            temp->setLabel(sharedorder->getLabel());
            temp->setGroup(groups[i]);
            indexLookup[groups[i]] = i;
            lookup.push_back(temp);
        }
        
        map<string, int> groupComboToColumn = getGroupComb(groups); //makes  'uniqueAB         uniqueAC  uniqueBC' if your groups are A, B, C
        
        SharedCollectorsCurveData ccd; ccd.registerDisplays(displays); //adds a displays to ccd management

        //convert freq percentage to number
        int increment = 1;
        if (percentFreq < 1.0) {  increment = numSeqs * percentFreq;  }
        else { increment = percentFreq;  }
        
        for(int i=0;i<numSeqs;i++){
            
            if (m->getControl_pressed()) { break;  }
            
            //get first sample
            individual chosen = sharedorder->get(i);
            lookup[indexLookup[chosen.group]]->increment(chosen.binNumber);
            
            //calculate at 0 and the given increment
            if((i == 0) || (i+1) % increment == 0){
                ccd.updateSharedData(lookup, i+1, groupComboToColumn);
            }
        }
        
        if (m->getControl_pressed()) { for (int j = 0; j < lookup.size(); j++) {  delete lookup[j]; }   return 1;  }
        
        //calculate last label if you haven't already
        if(numSeqs % increment != 0){ ccd.updateSharedData(lookup, numSeqs, groupComboToColumn); }
        
        //resets output files
        for(int i=0;i<displays.size();i++){ displays[i]->reset(); }
        
        //memory cleanup
        for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "Collect", "getSharedCurve");
        exit(1);
    }
}
/**************************************************************************************/
map<string, int> Collect::getGroupComb(vector<string> mGroups) {
	string group;
	numGroupComb = 0;
    map<string, int> groupComboToColumn;
    
    int numGroups = mGroups.size();
	for (int i = 0; i < (numGroups - 1); i++) {
		for (int l = i+1; l < numGroups; l++) {
			group = mGroups[i] +"_"+ mGroups[l];
			groupComb.push_back(group);
            groupComboToColumn[group] = numGroupComb;
			numGroupComb++;
		}
	}
    
    for(int i=0;i<displays.size();i++){
        bool hasLciHci = displays[i]->hasLciHci();
        groupLabel = "";
        for (int s = 0; s < groupComb.size(); s++) {
            if (hasLciHci) {  groupLabel +=  label +"_"+ groupComb[s] + "\t" + label + groupComb[s] + "lci\t" + label + groupComb[s] + "hci\t"; }
            else{  groupLabel += label +"_"+ groupComb[s] + "\t";  }
        }
        
        string groupLabelAll = groupLabel + label +"_"+ "all\t";
        if ((displays[i]->isCalcMultiple() ) && (displays[i]->getAll() )) { displays[i]->init(groupLabelAll); }
        else {  displays[i]->init(groupLabel);  }
    }
    
    return groupComboToColumn;
}
/**************************************************************************************/
