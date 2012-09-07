//
//  subsample.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "subsample.h"
//**********************************************************************************************************************
Tree* SubSample::getSample(Tree* T, CountTable* ct, CountTable* newCt, int size) {
    try {
        Tree* newTree = NULL;
        
        //remove seqs not in sample from counttable
        vector<string> Groups = ct->getNamesOfGroups();
        newCt->copy(ct); 
        newCt->addGroup("doNotIncludeMe");
        
        map<string, int> doNotIncludeTotals; 
        vector<string> namesSeqs = ct->getNamesOfSeqs();
        for (int i = 0; i < namesSeqs.size(); i++) {  doNotIncludeTotals[namesSeqs[i]] = 0; }
    
        for (int i = 0; i < Groups.size(); i++) {
            if (m->inUsersGroups(Groups[i], m->getGroups())) {
                if (m->control_pressed) { break; }
        
                int thisSize = ct->getGroupCount(Groups[i]);
                
                if (thisSize >= size) {	
                    
                    vector<string> names = ct->getNamesOfSeqs(Groups[i]);
                    vector<int> random;
                    for (int j = 0; j < names.size(); j++) {
                        int num = ct->getGroupCount(names[j], Groups[i]);
                        for (int k = 0; k < num; k++) { random.push_back(j); }
                    }
                    random_shuffle(random.begin(), random.end());
                    
                    vector<int> sampleRandoms; sampleRandoms.resize(names.size(), 0);
                    for (int j = 0; j < size; j++) { sampleRandoms[random[j]]++; }
                    for (int j = 0; j < sampleRandoms.size(); j++) {
                        newCt->setAbund(names[j], Groups[i], sampleRandoms[j]);
                    }
                    sampleRandoms.clear(); sampleRandoms.resize(names.size(), 0);
                    for (int j = size; j < thisSize; j++) { sampleRandoms[random[j]]++; }
                    for (int j = 0; j < sampleRandoms.size(); j++) {  doNotIncludeTotals[names[j]] += sampleRandoms[j]; }
                }else {  m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->control_pressed = true; }
            }

        }
        
        for (map<string, int>::iterator it = doNotIncludeTotals.begin(); it != doNotIncludeTotals.end(); it++) {  
            newCt->setAbund(it->first, "doNotIncludeMe", it->second);
        } 
        
        newTree = new Tree(newCt);
        newTree->getCopy(T, true);
        
        return newTree;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getSample-Tree");
        exit(1);
    }
}
//**********************************************************************************************************************
//assumes whole maps dupName -> uniqueName
map<string, string> SubSample::deconvolute(map<string, string> whole, vector<string>& wanted) {
    try {
        map<string, string> nameMap;
        
        //whole will be empty if user gave no name file, so we don't need to make a new one
        if (whole.size() == 0) { return nameMap; }
        
        vector<string> newWanted;
        for (int i = 0; i < wanted.size(); i++) {
            
            if (m->control_pressed) { break; }
            
            string dupName = wanted[i];
            
            map<string, string>::iterator itWhole = whole.find(dupName);
            if (itWhole != whole.end()) {
                string repName = itWhole->second;
                
                //do we already have this rep?
                map<string, string>::iterator itName = nameMap.find(repName);
                if (itName != nameMap.end()) { //add this seqs to dups list
                    (itName->second) += "," + dupName;
                }else { //first sighting of this seq
                    nameMap[repName] = dupName;
                    newWanted.push_back(repName);
                }
            }else { m->mothurOut("[ERROR]: "+dupName+" is not in your name file, please correct.\n"); m->control_pressed = true; }
        }
        
        wanted = newWanted;
        return nameMap;
    }
	catch(exception& e) {
		m->errorOut(e, "SubSample", "deconvolute");
		exit(1);
	}
}
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
		m->errorOut(e, "SubSample", "getSample-shared");
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
int SubSample::getSample(SAbundVector*& sabund, int size) {
	try {
		
        OrderVector* order = new OrderVector();
        *order = sabund->getOrderVector(NULL);
        
		int numBins = order->getNumBins();
		int thisSize = order->getNumSeqs();
        
		if (thisSize > size) {
			random_shuffle(order->begin(), order->end());
			
            RAbundVector* rabund = new RAbundVector(numBins);
			rabund->setLabel(order->getLabel());

			for (int j = 0; j < size; j++) {
                
				if (m->control_pressed) { delete order; delete rabund; return 0; }
				
				int bin = order->get(j);
				
				int abund = rabund->get(bin);
				rabund->set(bin, (abund+1));
			}
			
            delete sabund;
            sabund = new SAbundVector();
            *sabund = rabund->getSAbundVector();
            delete rabund;
            
		}else if (thisSize < size) { m->mothurOut("[ERROR]: The size you requested is larger than the number of sequences in the sabund vector. You requested " + toString(size) + " and you only have " + toString(thisSize) + " seqs in your sabund vector.\n"); m->control_pressed = true; }
		
		if (m->control_pressed) { return 0; }
        
		delete order;
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSample");
		exit(1);
	}
}			
//**********************************************************************************************************************


