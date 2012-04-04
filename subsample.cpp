//
//  subsample.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "subsample.h"

//**********************************************************************************************************************
vector<string> SubSample::getSample(vector<SharedRAbundVector*>& thislookup, int size) {
	try {
		
		//save mothurOut's binLabels to restore for next label
		vector<string> saveBinLabels = m->currentBinLabels;
		
		int numBins = thislookup[0]->getNumBins();
		for (int i = 0; i < thislookup.size(); i++) {		
			int thisSize = thislookup[i]->getNumSeqs();
			
			if (thisSize != size) {
				
				string thisgroup = thislookup[i]->getGroup();
				
				OrderVector order;
				for(int p=0;p<numBins;p++){
					for(int j=0;j<thislookup[i]->getAbundance(p);j++){
						order.push_back(p);
					}
				}
				random_shuffle(order.begin(), order.end());
				
				SharedRAbundVector* temp = new SharedRAbundVector(numBins);
				temp->setLabel(thislookup[i]->getLabel());
				temp->setGroup(thislookup[i]->getGroup());
				
				delete thislookup[i];
				thislookup[i] = temp;
				
				
				for (int j = 0; j < size; j++) {
					
					if (m->control_pressed) {  return m->currentBinLabels; }
					
					int bin = order.get(j);
					
					int abund = thislookup[i]->getAbundance(bin);
					thislookup[i]->set(bin, (abund+1), thisgroup);
				}	
			}
		}
		
		//subsampling may have created some otus with no sequences in them
		eliminateZeroOTUS(thislookup);
		
		if (m->control_pressed) { return m->currentBinLabels; }
		
		//save mothurOut's binLabels to restore for next label
        vector<string> subsampleBinLabels = m->currentBinLabels;
		m->currentBinLabels = saveBinLabels;
		
		return subsampleBinLabels;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSample", "getSample");
		exit(1);
	}
}	
//**********************************************************************************************************************
int SubSample::eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup) {
	try {
		
		vector<SharedRAbundVector*> newLookup;
		for (int i = 0; i < thislookup.size(); i++) {
			SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(thislookup[i]->getLabel());
			temp->setGroup(thislookup[i]->getGroup());
			newLookup.push_back(temp);
		}
		
		//for each bin
		vector<string> newBinLabels;
		string snumBins = toString(thislookup[0]->getNumBins());
		for (int i = 0; i < thislookup[0]->getNumBins(); i++) {
			if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
			
			//look at each sharedRabund and make sure they are not all zero
			bool allZero = true;
			for (int j = 0; j < thislookup.size(); j++) {
				if (thislookup[j]->getAbundance(i) != 0) { allZero = false;  break;  }
			}
			
			//if they are not all zero add this bin
			if (!allZero) {
				for (int j = 0; j < thislookup.size(); j++) {
					newLookup[j]->push_back(thislookup[j]->getAbundance(i), thislookup[j]->getGroup());
				}
				//if there is a bin label use it otherwise make one
				string binLabel = "Otu";
				string sbinNumber = toString(i+1);
				if (sbinNumber.length() < snumBins.length()) { 
					int diff = snumBins.length() - sbinNumber.length();
					for (int h = 0; h < diff; h++) { binLabel += "0"; }
				}
				binLabel += sbinNumber; 
				if (i < m->currentBinLabels.size()) {  binLabel = m->currentBinLabels[i]; }
				
				newBinLabels.push_back(binLabel);
			}
		}
		
		for (int j = 0; j < thislookup.size(); j++) {  delete thislookup[j];  }
		thislookup.clear();
		
		thislookup = newLookup;
		m->currentBinLabels = newBinLabels;
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSample", "eliminateZeroOTUS");
		exit(1);
	}
}


//**********************************************************************************************************************


